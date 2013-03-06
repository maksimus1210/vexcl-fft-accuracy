#include <cmath>
#include <cassert>
#include <array>
#include <complex>
#include <limits>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>

int power_of_two(unsigned int n) {
    return (((n) > 0) && (((n) & ((n) - 1)) == 0));
}

static int pow2_atleast(int x) {
    int h;
    for(h = 1; h < x; h = 2 * h);
    return h;
}



template<class N = boost::multiprecision::static_mpfr_float_50>
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
                const N alpha = pi * N(2) * N(j) / N(n);
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
        for(int i = 0; i < n; ++i) b[i] = std::conj(w[i]) * a[i];
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

    N compare(std::vector<CN> &a, std::vector<CN> &b, double norm) {
        const bool inf = norm == std::numeric_limits<double>::infinity();
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

    template<class T>
    double accuracy(const std::vector<std::complex<T>> &a, const std::vector<std::complex<T>> &ffta, int sign, bool forward = true, double norm = 2) {
        const int n = a.size();
        assert(n == ffta.size());
        std::vector<CN> b(n), fftb(n);
        for(int i = 0 ; i < n ; i++) {
            b[i] = CN(a[i]);
            fftb[i] = CN(ffta[i]);
        }
        if(forward) { // forward error
            fft1(b, sign);
        } else { // backward error
            for(int i = 0; i < n; ++i) fftb[i] /= n;
            fft1(fftb, -sign);
        }
        return static_cast<double>(compare(b, fftb, norm));
    }

};
