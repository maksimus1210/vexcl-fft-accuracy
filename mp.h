#include <cmath>
#include <cassert>
#include <array>
#include <complex>
#include <limits>
#include <vector>

#include <boost/math/constants/constants.hpp>

int power_of_two(unsigned int n) {
    return (((n) > 0) && (((n) & ((n) - 1)) == 0));
}

static int pow2_atleast(int x) {
    int h;
    for(h = 1; h < x; h = 2 * h);
    return h;
}

template<class N>
struct simple_fft {
    typedef std::complex<N> CN;

    std::vector<CN> w, y;
    simple_fft() {}

    /// exp(i 2 pi m / n)
    static CN exp(int m, int n) {
        static int cached_n = -1;
        static CN w[64];
        if(n != cached_n) {
            const N pi = boost::math::constants::pi<N>();
            for(int j = 1, k = 0; j < n; j += j, ++k) {
                const N alpha = pi * 2 * j / n;
                w[k] = CN(cos(alpha), sin(alpha));
            }
            cached_n = n;
        }
        CN v(1.0, 0.0);
        if(m > 0) {
            for(int k = 0; m; ++k, m >>= 1)
                if(m & 1) v *= w[k];
        } else {
            m = -m;
            for(int k = 0; m; ++k, m >>= 1)
                if(m & 1) v *= std::conj(w[k]);
        }
        return v;
    }

    void bitrev(std::vector<CN> &a) {
        const int n = a.size();
        for(int i = 0, j = 0; i < n - 1; ++i) {
            if(i < j) std::swap(a[i], a[j]);
            // bit reversed counter
            int m = n;
            do {
                m >>= 1;
                j ^= m;
            } while(!(j & m));
        }
    }

    void fft0(std::vector<CN> &a, int sign) {
        const int n = a.size();
        bitrev(a);
        for(int i = 1; i < n; i = 2 * i) {
            for(int j = 0; j < i; ++j) {
                const CN w = exp(sign * j, 2 * i);
                for(int k = j; k < n; k += 2 * i) {
                    const CN v = a[k];
                    const CN x = a[k + i] * w;
                    a[k] = v + x;
                    a[k + i] = v - x;
                }
            }
        }
    }

    void bluestein(std::vector<CN> &a) {
        const int n = a.size();
        const int nb = pow2_atleast(2 * n);
        std::vector<CN> b(nb);
        if(w.size() != n) {
            w.resize(n);
            // w[k] = exp(i 2 pi k^2 / (2 n))
            int ksq = 1; // (-1)^2
            for (int k = 0; k < n; ++k) {
                // careful with overflow
                ksq = ksq + 2 * k - 1;
                while(ksq > 2 * n) ksq -= 2 * n;
                w[k] = exp(ksq, 2 * n);
            }
            y.resize(nb);
            for(int i = 0; i < nb; ++i) y[i] = CN(0.0, 0.0);
            for(int i = 0; i < n; ++i) y[i] = w[i];
            for(int i = 1; i < n; ++i) y[nb-i] = w[i];
            fft0(y, -1);
        }
        for(int i = 0; i < nb; ++i) b[i] = CN(0.0, 0.0);
        for(int i = 0; i < n; ++i) b[i] = conj(w[i]) * a[i];
        // scaled convolution b * y
        fft0(b, -1);
        for(int i = 0; i < nb; ++i) b[i] *= y[i];
        fft0(b, 1);
        for(int i = 0; i < n; ++i) a[i] = conj(w[i]) * b[i] / N(nb);
    }

    void swapri(std::vector<CN> &a) {
        for(auto &x : a) x = CN(std::imag(x), std::real(x));
    }

    void fft1(std::vector<CN> &a, int sign) {
        if(power_of_two(a.size()))
            fft0(a, sign);
        else {
            if(sign == 1) swapri(a);
            bluestein(a);
            if(sign == 1) swapri(a);
        }
    }

    template<class T>
    std::vector<CN> operator()(const std::vector<T> &a, int sign = -1) {
        const int n = a.size();
        std::vector<CN> b(a.begin(), a.end());
        fft1(b, sign);
        if(sign == 1) for(auto &x : b) x /= n;
        return b;
    }

    template<class T1, class T2>
    N compare(const std::vector<T1> &_a, const std::vector<T2> &_b, double norm = 2) const {
        const bool inf = norm == std::numeric_limits<double>::infinity();
        std::vector<CN> a(_a.begin(), _a.end()), b(_b.begin(), _b.end());
        N e(0), m(0);
        for(int i = 0; i < a.size(); ++i) {
#define NORM(b, a) {N __v = abs(a); b = inf ? std::max(b, __v) : (b + pow(__v, norm));}
            NORM(m, std::real(a[i]));
            NORM(m, std::imag(a[i]));
            const CN d = a[i] - b[i];
            NORM(e, std::real(d));
            NORM(e, std::imag(d));
#undef NORM
        }
        if(inf) return e / m;
        else return pow(e / m, 1 / norm);
    }
};
