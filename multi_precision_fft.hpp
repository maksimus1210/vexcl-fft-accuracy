#ifndef MULTI_PRECISION_FFT_HPP
#define MULTI_PRECISION_FFT_HPP

/* Based on FFTw's benchfft-3.1/libbench/mp.c */

#include <complex>
#include <vector>
#include <boost/math/constants/constants.hpp>

inline bool is_pow_2(size_t n) {
    return n & (n - 1) == 0;
}

inline size_t next_pow_2(size_t x) {
    size_t h = 1;
    while(h < x) h *= 2;
    return h;
}

/// exp(i 2 pi m / n)
template<class N>
static std::complex<N> exp(int m, int n) {
    static const N pi = boost::math::constants::pi<N>();
    return std::polar<N>(1, pi * 2 * m / n);
}

/// bit-reverse transpose the array
template<class T>
static void bitrev(std::vector<T> &a) {
    const size_t n = a.size();
    for(size_t i = 0, j = 0 ; i < n - 1 ; i++) {
        if(i < j) std::swap(a[i], a[j]);
        // bit reversed counter
        size_t m = n;
        do {
            m >>= 1;
            j ^= m;
        } while(!(j & m));
    }
}

/// compute a power-of-2 FFT in-place.
template<class N>
static void fft0(std::vector<std::complex<N>> &a, int sign) {
    const size_t n = a.size();
    bitrev(a);
    for(size_t i = 1 ; i < n ; i *= 2)
        for(size_t j = 0 ; j < i ; j++) {
            const auto w = exp<N>(sign * j, 2 * i);
            for(size_t k = j ; k < n ; k += 2 * i) {
                const auto v = a[k];
                const auto x = a[k + i] * w;
                a[k] = v + x;
                a[k + i] = v - x;
            }
        }
}

/// compute a Bluestein Chirp-Z transform in-place.
template<class N>
void bluestein(std::vector<std::complex<N>> &a) {
    const size_t n = a.size();
    const size_t nb = next_pow_2(2 * n);
    const std::complex<N> zero(0, 0);
    std::vector<std::complex<N>> w(n), y(nb, zero), b(nb, zero);
    for(size_t k = 0 ; k < n ; k++)
        w[k] = exp<N>(((unsigned long long)k * k) % (2 * n), 2 * n);
    for(size_t i = 0 ; i < n ; i++) y[i] = w[i];
    for(size_t i = 1 ; i < n ; i++) y[nb - i] = w[i];
    fft0(y, -1);
    for(size_t i = 0 ; i < n ; i++) b[i] = conj(w[i]) * a[i];
    // scaled convolution b * y
    fft0(b, -1);
    for(size_t i = 0 ; i < nb ; i++) b[i] *= y[i];
    fft0(b, 1);
    for(size_t i = 0 ; i < n ; i++) a[i] = conj(w[i]) * b[i] / N(nb);
}

/// Swap real and imaginary parts of all elements.
template<class N>
static void swapri(std::vector<std::complex<N>> &a) {
    for(auto &x : a) x = std::complex<N>(std::imag(x), std::real(x));
}

/// Compute an arbitrary-sized FFT in-place.
template<class N>
static void fft1(std::vector<std::complex<N>> &a, int sign) {
    if(is_pow_2(a.size()))
        fft0(a, sign);
    else {
        if(sign == 1) swapri(a);
        bluestein(a);
        if(sign == 1) swapri(a);
    }
}

/// Compute an arbitrary-sized FFT out-of-place.
template<class N, class T>
static std::vector<std::complex<N>> simple_fft(const std::vector<T> &a, int sign = -1) {
    const int n = a.size();
    std::vector<std::complex<N>> b(a.begin(), a.end());
    fft1(b, sign);
    if(sign == 1) for(auto &x : b) x /= n;
    return b;
}

/// Compute the difference between two arrays using a given norm, precision determined by first argument.
template<class T1, class T2>
static T1 compare(const std::vector<std::complex<T1>> &_a, const std::vector<std::complex<T2>> &_b, double norm = 2) {
    const bool inf = norm == std::numeric_limits<double>::infinity();
    std::vector<std::complex<T1>> a(_a.begin(), _a.end()), b(_b.begin(), _b.end());
    T1 e(0), m(0);
    for(int i = 0; i < a.size(); ++i) {
#define NORM(b, a) {T1 __v = abs(a); b = inf ? std::max(b, __v) : (b + pow(__v, norm));}
        NORM(m, std::real(a[i]));
        NORM(m, std::imag(a[i]));
        const auto d = a[i] - b[i];
        NORM(e, std::real(d));
        NORM(e, std::imag(d));
#undef NORM
    }
    if(inf) return e / m;
    else return pow(e / m, 1 / norm);
}


#endif
