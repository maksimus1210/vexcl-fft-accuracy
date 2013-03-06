#include <cmath>
#include <cassert>
#include <array>
#include <complex>
#include <limits>
#include <vector>

typedef unsigned short DG;
typedef unsigned long ACC;
const size_t BITS_IN_REAL = 53; // mantissa

const size_t SHFT = 16;
const size_t RADIX = 65536L;
#define IRADIX (1.0 / RADIX)
#define LO(x) ((x) & (RADIX - 1))
#define HI(x) ((x) >> SHFT)
#define HI_SIGNED(x) \
   ((((x) + (ACC)(RADIX >> 1) * RADIX) >> SHFT) - (RADIX >> 1))
const short ZEROEXP = -32768;

const size_t LEN = 10;

struct N {
    short sign, expt;
    std::array<DG, LEN> d;
    N(short s, short e, std::array<DG,LEN> dd) : sign(s), expt(e), d(dd) {}

    N(int x = 0) : sign(1), expt(1) {
        for(size_t i = 0 ; i < LEN ; i++) d[i] = 0;
        if(x == 0) {
            expt = ZEROEXP;
        } else {
            if(x < 0) {
                x = -x;
                sign = -1;
            }
            d[LEN - 1] = x;
        }
    }

    template<class T>
    N(T x) {
        *this = N();
        if(x == 0.0) return;
        if(x < 0) { sign = -1; x = -x; }
        expt = 0;
        while(x >= 1.0) { x *= IRADIX; expt++; }
        while(x < IRADIX) { x *= RADIX; expt--; }
        for(int i = LEN - 1; i >= 0 && x != 0.0; --i) {
            x *= RADIX;
            T y = floor(x);
            d[i] = (DG)y;
            x -= y;
        }
    }

    N(DG *_d, int e, int s, int l) {
        *this = N();
        int i;
        for(i = l - 1; i >= 0; --i, --e)
            if(_d[i] != 0) break;
        if(i < 0) return;
        expt = e;
        sign = s;
        if(i >= LEN - 1) {
            for(int j = LEN - 1; j >= 0; --i, --j) d[j] = _d[i];
        } else {
            int j;
            for(j = LEN - 1; i >= 0; --i, --j) d[j] = _d[i];
            for(; j >= 0; --j) d[j] = 0;
        }
    }
};


/// compare absolute values
static int abscmp(const N &a, const N &b) {
    if(a.expt > b.expt) return 1;
    if(a.expt < b.expt) return -1;
    for(int i = LEN - 1; i >= 0; --i) {
        if(a.d[i] > b.d[i]) return 1;
        if(a.d[i] < b.d[i]) return -1;
    }
    return 0;
}

static bool operator==(const N &a, const N &b) {
    return (a.sign == b.sign) && (abscmp(a, b) == 0);
}

/// add magnitudes, for |a| >= |b|
static N addmag0(int s, const N &a, const N &b) {
    int ia, ib;
    ACC r = 0;
    DG d[LEN + 1];
    for(ia = 0, ib = a.expt - b.expt; ib < LEN; ++ia, ++ib) {
        r += (ACC)a.d[ia] + (ACC)b.d[ib];
        d[ia] = LO(r);
        r = HI(r);
    }
    for(; ia < LEN; ++ia) {
        r += (ACC)a.d[ia];
        d[ia] = LO(r);
        r = HI(r);
    }
    d[ia] = LO(r);
    return N(d, a.expt + 1, s * a.sign, LEN + 1);
}

static N addmag(int s, const N &a, const N &b) {
    return (abscmp(a, b) > 0) ? addmag0(1, a, b) : addmag0(s, b, a);
}

/// subtract magnitudes, for |a| >= |b|
static N submag0(int s, const N &a, const N &b) {
    int ia, ib;
    ACC r = 0;
    DG d[LEN];
    for (ia = 0, ib = a.expt - b.expt; ib < LEN; ++ia, ++ib) {
        r += (ACC)a.d[ia] - (ACC)b.d[ib];
        d[ia] = LO(r);
        r = HI_SIGNED(r);
    }
    for (; ia < LEN; ++ia) {
        r += (ACC)a.d[ia];
        d[ia] = LO(r);
        r = HI_SIGNED(r);
    }
    return N(d, a.expt, s * a.sign, LEN);
}

static N submag(int s, const N a, const N &b) {
    return (abscmp(a, b) > 0) ? submag0(1, a, b) : submag0(s, b, a);
}

static N operator+(const N &a, const N &b) {
    return (a.sign == b.sign) ? addmag(1, a, b) : submag(1, a, b);
}

static N operator+=(N &a, const N &b) {
    a = a + b;
    return a;
}

static N operator-(const N &a, const N &b) {
    return (a.sign == b.sign) ? submag(-1, a, b) : addmag(-1, a, b);
}

static N operator-=(N &a, const N &b) {
    a = a - b;
    return a;
}

static N operator*(const N &a, const N &b) {
    DG d[2 * LEN];
    int i, j, k;
    ACC r;
    for(i = 0; i < LEN; ++i)
        d[2 * i] = d[2 * i + 1] = 0;
    for(i = 0; i < LEN; ++i) {
        ACC ai = a.d[i];
        if(ai) {
            r = 0;
            for(j = 0, k = i; j < LEN; ++j, ++k) {
                r += ai * (ACC)b.d[j] + (ACC)d[k];
                d[k] = LO(r);
                r = HI(r);
            }
            d[k] = LO(r);
        }
    }
    return N(d, a.expt + b.expt, a.sign * b.sign, 2 * LEN);
}

static N operator*=(N &a, const N &b) {
    a = a * b;
    return a;
}

template<typename T>
static T to(const N &a) {
    T h, l, f;
    int i, bits;
    ACC r;
    DG sticky;

    if (a.expt != ZEROEXP) {
        f = IRADIX;
        i = LEN;

        bits = 0;
        h = (r = a.d[--i]) * f; f *= IRADIX;
        for (bits = 0; r > 0; ++bits)
            r >>= 1;

        // first digit
        while (bits + SHFT <= BITS_IN_REAL) {
            h += a.d[--i] * f;
            f *= IRADIX;
            bits += SHFT;
        }

        // guard digit (leave one bit for sticky bit, hence `<' instead of `<=')
        bits = 0; l = 0.0;
        while (bits + SHFT < BITS_IN_REAL) {
            l += a.d[--i] * f;
            f *= IRADIX;
            bits += SHFT;
        }

        // sticky bit
        sticky = 0;
        while(i > 0) sticky |= a.d[--i];
        if(sticky) l += (RADIX / 2) * f;

        h += l;

        for(i = 0; i < a.expt; ++i) h *= (T)RADIX;
        for(i = 0; i > a.expt; --i) h *= IRADIX;
        if(a.sign == -1) h = -h;
        return h;
    } else {
        return 0.0;
    }
}

static N operator-(const N &a) {
    N b = a;
    b.sign = -b.sign;
    return b;
}

static N inv(const N &a) {
    static const N one(1), two(2);
    N x(1.0 / to<double>(a)); // initial guess
    for(;;) { // Newton
        N z = two - a * x;
        if(one == z) break;
        x *= z;
    }
    return x;
}


static N sin(const N &a) {
    // 1 / 31!
    static const N i31fac(1, -7, {28087, 45433, 51357, 24545, 14291, 3954, 57879, 8109, 38716, 41382});
    N g = i31fac, b = g;
    const N a2 = a * a;
    for(int i = 31; i > 1; i -= 2) { // Taylor
        g *= N(i * (i - 1));
        b = g - a2 * b;
    }
    return a * b;
}

static N cos(const N &a) {
    // 1 / 32!
    static const N i32fac(1, -7, {52078, 60811, 3652, 39679, 37310, 47227, 28432, 57597, 13497, 1293});
    N g = i32fac, b = g;
    const N a2 = a * a;
    for(int i = 32; i > 0; i -= 2) { // Taylor
        g *= N(i * (i - 1));
        b = g - a2 * b;
    }
    return b;
}

/// 2 pi m / n
static N by2pi(double m, double n) {
    static const N n2pi(1,  1, {18450, 59017, 1760, 5212, 9779, 4518, 2886, 54545, 18558, 6});
    return N(m) * inv(N(n)) * n2pi;
}

static N sin2pi(double m, double n);

static N cos2pi(double m, double n) {
    if(m < 0) return cos2pi(-m, n);
    else if(m > n * 0.5) return cos2pi(n - m, n);
    else if(m > n * 0.25) return -sin2pi(m - n * 0.25, n);
    else if(m > n * 0.125) return sin2pi(n * 0.25 - m, n);
    else return cos(by2pi(m, n));
}

static N sin2pi(double m, double n) {
    if(m < 0) return -sin2pi(-m, n);
    else if(m > n * 0.5) return -sin2pi(n - m, n);
    else if(m > n * 0.25) return cos2pi(m - n * 0.25, n);
    else if(m > n * 0.125) return cos2pi(n * 0.25 - m, n);
    else return sin(by2pi(m, n));
}

// FFT stuff

struct CN : std::complex<N> {
    // from floating point
    template<class T>
    CN(const std::complex<T> &v) : std::complex<N>(N(std::real(v)), N(std::imag(v))) {}

    // polar form exp(i 2 pi m / n)
    CN(int m, int n) : std::complex<N>(cos2pi(m,n), sin2pi(m,n)) {}

    CN(const N &r = N(), const N &i = N()) : std::complex<N>(r, i) {}
};

int power_of_two(unsigned int n) {
    return (((n) > 0) && (((n) & ((n) - 1)) == 0));
}

static int pow2_atleast(int x) {
    int h;
    for(h = 1; h < x; h = 2 * h);
    return h;
}

/// exp(i 2 pi m / n)
static std::complex<N> exp(int m, int n) {
    static int cached_n = -1;
    static std::complex<N> w[64];
    if(n != cached_n) {
        for(int j = 1, k = 0; j < n; j += j, ++k)
            w[k] = CN(j, n);
        cached_n = n;
    }
    CN v(N(1));
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

struct simple_fft {
    std::vector<CN> w, y;
    simple_fft() {}

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
        const N nbinv = N(1.0 / nb); // exact because nb = 2^k
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
            for(int i = 0; i < nb; ++i) y[i] = CN();
            for(int i = 0; i < n; ++i) y[i] = w[i];
            for(int i = 1; i < n; ++i) y[nb-i] = w[i];
            fft0(y, -1);
        }
        for(int i = 0; i < nb; ++i) b[i] = CN();
        for(int i = 0; i < n; ++i) b[i] = std::conj(w[i]) * a[i];
        // scaled convolution b * y
        fft0(b, -1);
        for(int i = 0; i < nb; ++i) b[i] *= y[i];
        fft0(b, 1);
        for(int i = 0; i < n; ++i) a[i] = conj(w[i]) * b[i] * nbinv;
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

    double compare(std::vector<CN> &a, std::vector<CN> &b, double norm) {
        const double inf = norm == std::numeric_limits<double>::infinity();
        double e = 0, m = 0;
        for(int i = 0; i < a.size(); ++i) {
#define NORM(b, a) \
    {double __v = std::abs(to<double>(a)); \
    if(inf) b = std::max(b, __v); \
    else b += std::pow(__v, norm);}
            NORM(m, std::real(a[i]));
            NORM(m, std::imag(a[i]));
            const CN d = a[i] - b[i];
            NORM(e, std::real(d));
            NORM(e, std::imag(d));
#undef NORM
        }
        if(inf) return e / m;
        else return std::pow(e / m, 1 / norm);
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
            const N ninv = inv(N(n));
            for(int i = 0; i < n; ++i) fftb[i] *= ninv;
            fft1(fftb, -sign);
        }
        return compare(b, fftb, norm);
    }

};
