#include <cmath>
#include <cstdlib>

typedef T bench_real;
typedef T2 bench_complex;

int power_of_two(unsigned int n)
{
        return (((n) > 0) && (((n) & ((n) - 1)) == 0));
}

#define c_re(x) (x).s[0]
#define c_im(x) (x).s[1]

#define bench_malloc malloc
#define bench_free free

#define DG unsigned short
#define ACC unsigned long
#define REAL bench_real
#define BITS_IN_REAL 53 /* mantissa */

#define SHFT 16
#define RADIX 65536L
#define IRADIX (1.0 / RADIX)
#define LO(x) ((x) & (RADIX - 1))
#define HI(x) ((x) >> SHFT)
#define HI_SIGNED(x) \
   ((((x) + (ACC)(RADIX >> 1) * RADIX) >> SHFT) - (RADIX >> 1))
#define ZEROEXP (-32768)

#define LEN 10

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

    N(REAL x) {
        *this = N();
        if(x == 0.0) return;
        if(x < 0) { sign = -1; x = -x; }
        expt = 0;
        while(x >= 1.0) { x *= IRADIX; expt++; }
        while(x < IRADIX) { x *= RADIX; expt--; }
        for(int i = LEN - 1; i >= 0 && x != 0.0; --i) {
            x *= RADIX;
            REAL y = floor(x);
            d[i] = (DG)y;
            x -= y;
        }
    }
};


static N pack(DG *d, int e, int s, int l) {
    int i;
    for(i = l - 1; i >= 0; --i, --e)
        if(d[i] != 0) break;
    if(i < 0) return N(0);
    N a;
    a.expt = e;
    a.sign = s;
    if(i >= LEN - 1) {
        for(int j = LEN - 1; j >= 0; --i, --j) a.d[j] = d[i];
    } else {
        int j;
        for(j = LEN - 1; i >= 0; --i, --j) a.d[j] = d[i];
        for(; j >= 0; --j) a.d[j] = 0;
    }
    return a;
}


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
    return pack(d, a.expt + 1, s * a.sign, LEN + 1);
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
    return pack(d, a.expt, s * a.sign, LEN);
}

static N submag(int s, const N a, const N &b) {
    return (abscmp(a, b) > 0) ? submag0(1, a, b) : submag0(s, b, a);
}

static N operator+(const N &a, const N &b) {
    return (a.sign == b.sign) ? addmag(1, a, b) : submag(1, a, b);
}

static N operator-(const N &a, const N &b) {
    return (a.sign == b.sign) ? submag(-1, a, b) : addmag(-1, a, b);
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
    return pack(d, a.expt + b.expt, a.sign * b.sign, 2 * LEN);
}

static N operator*=(N &a, const N &b) {
    a = a * b;
    return a;
}

static REAL toreal(const N &a) {
    REAL h, l, f;
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

        for(i = 0; i < a.expt; ++i) h *= (REAL)RADIX;
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
    const N one(1), two(2);
    N x(1.0 / toreal(a)); // initial guess
    for(;;) { // Newton
        N z = two - a * x;
        if(one == z) break;
        x *= z;
    }
    return x;
}


// 2 pi
static const N n2pi(1,  1, {18450, 59017, 1760, 5212, 9779, 4518, 2886, 54545, 18558, 6});

// 1 / 31!
static const N i31fac(1, -7, {28087, 45433, 51357, 24545, 14291, 3954, 57879, 8109, 38716, 41382});

// 1 / 32!
static const N i32fac(1, -7, {52078, 60811, 3652, 39679, 37310, 47227, 28432, 57597, 13497, 1293});


static N sin(const N &a) {
    N g = i31fac, b = g;
    const N a2 = a * a;
    for(int i = 31; i > 1; i -= 2) { // Taylor
        g *= N(i * (i - 1));
        b = g - a2 * b;
    }
    return a * b;
}

static N cos(const N &a) {
    N g = i32fac, b = g;
    const N a2 = a * a;
    for(int i = 32; i > 0; i -= 2) { // Taylor
        g *= N(i * (i - 1));
        b = g - a2 * b;
    }
    return b;
}

/// 2 pi m / n
static N by2pi(REAL m, REAL n) {
    return N(m) * inv(N(n)) * n2pi;
}

static N sin2pi(REAL m, REAL n);
static N cos2pi(REAL m, REAL n) {
    if(m < 0) return cos2pi(-m, n);
    else if(m > n * 0.5) return cos2pi(n - m, n);
    else if(m > n * 0.25) return -sin2pi(m - n * 0.25, n);
    else if(m > n * 0.125) return sin2pi(n * 0.25 - m, n);
    else return cos(by2pi(m, n));
}

static N sin2pi(REAL m, REAL n) {
    if(m < 0) return -sin2pi(-m, n);
    else if(m > n * 0.5) return -sin2pi(n - m, n);
    else if(m > n * 0.25) return cos2pi(m - n * 0.25, n);
    else if(m > n * 0.125) return cos2pi(n * 0.25 - m, n);
    else return sin(by2pi(m, n));
}

// FFT stuff

struct CN {
    N r, i;
    CN(N r = N(), N i = N()) : r(r), i(i) {}

    CN operator*(const CN &v) const {
        return CN(r * v.r - i * v.i, r * v.i + i * v.r);
    }

    CN operator+(const CN &v) const {
        return CN(r + v.r, i + v.i);
    }

    CN operator-(const CN &v) const {
        return CN(r - v.r, i - v.i);
    }
};

static CN conj(const CN &v) {
    return CN(v.r, -v.i);
}

static CN exp(int m, int n) {
    static int cached_n = -1;
    static CN w[64];
    if(n != cached_n) {
        for(int j = 1, k = 0; j < n; j += j, ++k)
            w[k] = CN(cos2pi(j, n), sin2pi(j, n));
        cached_n = n;
    }
    CN v(N(1));
    if(m > 0) {
        for(int k = 0; m; ++k, m >>= 1)
            if(m & 1) v = w[k] * v;
    } else {
        m = -m;
        for(int k = 0; m; ++k, m >>= 1)
            if(m & 1) v = conj(w[k]) * v;
    }
    return v;
}

static void bitrev(int n, CN *a) {
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

static void fft0(int n, CN *a, int sign) {
    bitrev(n, a);
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

/// a[2*k]+i*a[2*k+1] = exp(2*pi*i*k^2/(2*n))
static void bluestein_sequence(int n, CN *a) {
    int ksq = 1; // (-1)^2
    for (int k = 0; k < n; ++k) {
        // careful with overflow
        ksq = ksq + 2 * k - 1;
        while(ksq > 2 * n) ksq -= 2 * n;
        a[k] = exp(ksq, 2 * n);
    }
}

static int pow2_atleast(int x) {
    int h;
    for(h = 1; h < x; h = 2 * h);
    return h;
}

static CN *cached_bluestein_w = 0;
static CN *cached_bluestein_y = 0;
static int cached_bluestein_n = -1;

static void bluestein(int n, CN *a) {
    const int nb = pow2_atleast(2 * n);
    CN *b = new CN[nb];
    CN *w = cached_bluestein_w;
    CN *y = cached_bluestein_y;
    const N nbinv = N(1.0 / nb); // exact because nb = 2^k
    if(cached_bluestein_n != n) {
        if(w) delete [] w;
        if(y) delete [] y;
        w = new CN[n];
        y = new CN[nb];
        bluestein_sequence(n, w);
        for(int i = 0; i < nb; ++i) y[i] = CN();
        for(int i = 0; i < n; ++i) y[i] = w[i];
        for(int i = 1; i < n; ++i) y[nb-i] = w[i];
        fft0(nb, y, -1);
        cached_bluestein_n = n;
        cached_bluestein_w = w;
        cached_bluestein_y = y;
    }
    for(int i = 0; i < nb; ++i) b[i] = CN(0);
    for(int i = 0; i < n; ++i) b[i] = conj(w[i]) * a[i];

    // scaled convolution b * y
    fft0(nb, b, -1);

    for(int i = 0; i < nb; ++i)
        b[i] = b[i] * y[i];
    fft0(nb, b, 1);

    for(int i = 0; i < n; ++i) {
        a[i] = conj(w[i]) * b[i];
        a[i].r *= nbinv;
        a[i].i *= nbinv;
    }
    delete [] b;
}

static void swapri(int n, CN *a) {
    for(int i = 0; i < n; ++i)
        std::swap(a[i].r, a[i].i);
}

static void fft1(int n, CN *a, int sign) {
    if(power_of_two(n))
        fft0(n, a, sign);
    else {
        if(sign == 1) swapri(n, a);
        bluestein(n, a);
        if(sign == 1) swapri(n, a);
    }
}

static void fromrealv(int n, const bench_complex *a, CN *b) {
    for (int i = 0; i < n; ++i) {
        b[i].r = N(c_re(a[i]));
        b[i].i = N(c_im(a[i]));
    }
}

static void compare(int n, const N *a, const N *b, double err[3]) {
    double e1 = 0, e2 = 0, einf = 0;
    double n1 = 0, n2 = 0, ninf = 0;
    for(int i = 0; i < 2 * n; ++i) {
        N dd = a[i] - b[i];
        {
            double d = std::abs(toreal(a[i]));
            n1 += d;
            n2 += d * d;
            if(d > ninf) ninf = d;
        }
        {
            double d = std::abs(toreal(dd));
            e1 += d;
            e2 += d * d;
            if(d > einf) einf = d;
        }
    }
    err[0] = e1 / n1;
    err[1] = sqrt(e2 / n2);
    err[2] = einf / ninf;
}

void fftaccuracy(int n, const bench_complex *a, const bench_complex *ffta, int sign, double err[6]) {
    CN *b = (CN *)bench_malloc(n * sizeof(CN));
    CN *fftb = (CN *)bench_malloc(n * sizeof(CN));
    N mn, ninv;
    int i;

    mn = N(n);
    ninv = inv(mn);

    // forward error
    fromrealv(n, a, b); fromrealv(n, ffta, fftb);
    fft1(n, b, sign);
    compare(n, (N*)b, (N*)fftb, err);

    // backward error
    fromrealv(n, a, b); fromrealv(n, ffta, fftb);
    for (i = 0; i < n; ++i) {fftb[i].r *= ninv; fftb[i].i *= ninv;}
    fft1(n, fftb, -sign);
    compare(n, (N*)b, (N*)fftb, err + 3);

    bench_free(fftb);
    bench_free(b);
}

void fftaccuracy_done(void)
{
    if (cached_bluestein_w) bench_free(cached_bluestein_w);
    if (cached_bluestein_y) bench_free(cached_bluestein_y);
}
