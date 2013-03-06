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
    N a,b;
    if (m < 0) a = cos2pi(-m, n);
    else if (m > n * 0.5) a = cos2pi(n - m, n);
    else if (m > n * 0.25) a = -sin2pi(m - n * 0.25, n);
    else if (m > n * 0.125) a = sin2pi(n * 0.25 - m, n);
    else { b = by2pi(m, n); a = cos(b); }
    return a;
}

static N sin2pi(REAL m, REAL n)
{
    N a,b;
    if (m < 0)  a = -sin2pi(-m, n);
    else if (m > n * 0.5) a = -sin2pi(n - m, n);
    else if (m > n * 0.25) a = cos2pi(m - n * 0.25, n);
    else if (m > n * 0.125) a = cos2pi(n * 0.25 - m, n);
    else {b = by2pi(m, n); a = sin(b);}
    return a;
}

/*----------------------------------------------------------------------*/
/* FFT stuff */

/* (r0 + i i0)(r1 + i i1) */
static void cmul(N &r0, N &i0, N &r1, N &i1, N &r2, N &i2)
{
    N s, t, q;
    s = r0 * r1;
    t = i0 * i1;
    q = s - t;
    s = r0 * i1;
    t = i0 * r1;
    i2 = s + t;
    r2 = q;
}

/* (r0 - i i0)(r1 + i i1) */
static void cmulj(N &r0, N &i0, N &r1, N &i1, N &r2, N &i2)
{
    N s, t, q;
    s = r0 * r1;
    t = i0 * i1;
    q = s + t;
    s = r0 * i1;
    t = i0 * r1;
    i2 = s - t;
    r2 = q;
}

static void cexp(int m, int n, N &r, N &i)
{
    static int cached_n = -1;
    static N w[64][2];
    int k, j;
    if (n != cached_n) {
        for (j = 1, k = 0; j < n; j += j, ++k) {
            w[k][0] = cos2pi(j, n);
            w[k][1] = sin2pi(j, n);
        }
        cached_n = n;
    }

    r = N(1);
    i = N(0);
    if (m > 0) {
        for (k = 0; m; ++k, m >>= 1) 
            if (m & 1)
                cmul(w[k][0], w[k][1], r, i, r, i);
    } else {
        m = -m;
        for (k = 0; m; ++k, m >>= 1) 
            if (m & 1)
                cmulj(w[k][0], w[k][1], r, i, r, i);
    }
}

static void bitrev(int n, N *a)
{
    int i, j, m;
    for (i = j = 0; i < n - 1; ++i) {
        if (i < j) {
            N t;
            t = a[2*i];
         a[2*i] = a[2*j];
         a[2*j] = t;
            t = a[2*i+1];
         a[2*i+1] = a[2*j+1];
         a[2*j+1] = t;
        }

        /* bit reversed counter */
        m = n; do { m >>= 1; j ^= m; } while (!(j & m));
    }
}

static void fft0(int n, N *a, int sign)
{
    int i, j, k;

    bitrev(n, a);
    for (i = 1; i < n; i = 2 * i) {
        for (j = 0; j < i; ++j) {
            N wr, wi;
            cexp(sign * (int)j, 2 * i, wr, wi);
            for (k = j; k < n; k += 2 * i) {
                N *a0 = a + 2 * k;
                N *a1 = a0 + 2 * i;
                N r0, i0, r1, i1, t0, t1, xr, xi;
                r0 = a0[0];
          i0 = a0[1];
                r1 = a1[0];
          i1 = a1[1];
                t0 = r1 * wr;
          t1 = i1 * wi;
          xr = t0 - t1;
                t0 = r1 * wi;
          t1 = i1 * wr;
          xi = t0 + t1;
                a0[0] = r0 + xr; 
          a0[1] = i0 + xi;
                a1[0] = r0 - xr;
          a1[1] = i0 - xi;
            }
        }
    }
}

/* a[2*k]+i*a[2*k+1] = exp(2*pi*i*k^2/(2*n)) */
static void bluestein_sequence(int n, N *a)
{
    int k, ksq, n2 = 2 * n;

    ksq = 1; /* (-1)^2 */
    for (k = 0; k < n; ++k) {
        /* careful with overflow */
        ksq = ksq + 2*k - 1; while (ksq > n2) ksq -= n2;
        cexp(ksq, n2, a[2*k], a[2*k+1]);
    }
}

static int pow2_atleast(int x)
{
    int h;
    for (h = 1; h < x; h = 2 * h)
        ;
    return h;
}

static N *cached_bluestein_w = 0;
static N *cached_bluestein_y = 0;
static int cached_bluestein_n = -1;

static void bluestein(int n, N *a)
{
    int nb = pow2_atleast(2 * n);
    N *b = (N *)bench_malloc(2 * nb * sizeof(N));
    N *w = cached_bluestein_w;
    N *y = cached_bluestein_y;
    N nbinv;
    int i;

    nbinv = N(1.0 / nb); /* exact because nb = 2^k */

    if (cached_bluestein_n != n) {
        if (w) bench_free(w);
        if (y) bench_free(y);
        w = (N *)bench_malloc(2 * n * sizeof(N));
        y = (N *)bench_malloc(2 * nb * sizeof(N));

        bluestein_sequence(n, w);
        for (i = 0; i < 2*nb; ++i)  y[i] = N(0);

        for (i = 0; i < n; ++i) {
            y[2*i] = w[2*i];
            y[2*i+1] = w[2*i+1];
        }
        for (i = 1; i < n; ++i) {
            y[2*(nb-i)] = w[2*i];
            y[2*(nb-i)+1] = w[2*i+1];
        }

        fft0(nb, y, -1);
        cached_bluestein_n = n;
        cached_bluestein_w = w;
        cached_bluestein_y = y;
    }

    for (i = 0; i < 2*nb; ++i)  b[i] = N(0);
    
    for (i = 0; i < n; ++i) 
        cmulj(w[2*i], w[2*i+1], a[2*i], a[2*i+1], b[2*i], b[2*i+1]);

    /* scaled convolution b * y */
    fft0(nb, b, -1);

    for (i = 0; i < nb; ++i) 
        cmul(b[2*i], b[2*i+1], y[2*i], y[2*i+1], b[2*i], b[2*i+1]);
    fft0(nb, b, 1);

    for (i = 0; i < n; ++i) {
        cmulj(w[2*i], w[2*i+1], b[2*i], b[2*i+1], a[2*i], a[2*i+1]);
        a[2*i] = nbinv * a[2*i];
        a[2*i+1] = nbinv * a[2*i+1];
    }

    bench_free(b);
}

static void swapri(int n, N *a)
{
    int i;
    for (i = 0; i < n; ++i) {
        N t;
        t = a[2 * i];
        a[2 * i] = a[2 * i + 1];
        a[2 * i + 1] = t;
    }
}

static void fft1(int n, N *a, int sign)
{
    if (power_of_two(n)) {
        fft0(n, a, sign);
    } else {
        if (sign == 1) swapri(n, a);
        bluestein(n, a);
        if (sign == 1) swapri(n, a);
    }
}

static void fromrealv(int n, bench_complex *a, N *b)
{
    int i;

    for (i = 0; i < n; ++i) {
        b[2 * i] = N(c_re(a[i]));
        b[2 * i + 1] = N(c_im(a[i]));
    }
}

static void compare(int n, N *a, N *b, double *err)
{
    int i;
    double e1, e2, einf;
    double n1, n2, ninf;

    e1 = e2 = einf = 0.0;
    n1 = n2 = ninf = 0.0;

#    define DO(x1, x2, xinf, var) {                   \
    double d = var;                              \
    if (d < 0) d = -d;                              \
    x1 += d; x2 += d * d; if (d > xinf) xinf = d;      \
}
        
    for (i = 0; i < 2 * n; ++i) {
        N dd;
        dd = a[i] - b[i];
        DO(n1, n2, ninf, toreal(a[i]));
        DO(e1, e2, einf, toreal(dd));
    }

#    undef DO
    err[0] = e1 / n1;
    err[1] = sqrt(e2 / n2);
    err[2] = einf / ninf;
}

void fftaccuracy(int n, bench_complex *a, bench_complex *ffta,
             int sign, double err[6])
{
    N *b = (N *)bench_malloc(2 * n * sizeof(N));
    N *fftb = (N *)bench_malloc(2 * n * sizeof(N));
    N mn, ninv;
    int i;

    mn = N(n);
    ninv = inv(mn);

    /* forward error */
    fromrealv(n, a, b); fromrealv(n, ffta, fftb);
    fft1(n, b, sign);
    compare(n, b, fftb, err);

    /* backward error */
    fromrealv(n, a, b); fromrealv(n, ffta, fftb);
    for (i = 0; i < 2 * n; ++i) fftb[i] = fftb[i] * ninv;
    fft1(n, fftb, -sign);
    compare(n, b, fftb, err + 3);

    bench_free(fftb);
    bench_free(b);
}

void fftaccuracy_done(void)
{
    if (cached_bluestein_w) bench_free(cached_bluestein_w);
    if (cached_bluestein_y) bench_free(cached_bluestein_y);
}
