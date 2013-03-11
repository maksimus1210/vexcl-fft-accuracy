#include "multi_precision_fft.hpp"

#include <random>
#include <boost/format.hpp>
#include <vexcl/vexcl.hpp>

#include <boost/multiprecision/mpfr.hpp>
typedef boost::multiprecision::static_mpfr_float_50 float50;
typedef boost::multiprecision::mpfr_float_1000 float1000;

#ifdef TEST_FFTW
#include <fftw3.h>
std::vector<std::complex<double>> fftw(const std::vector<std::complex<double>> &in) {
    std::vector<std::complex<double>> out(in.size());
    fftw_plan p = fftw_plan_dft_1d(in.size(),
        reinterpret_cast<fftw_complex *>(const_cast<std::complex<double> *>(in.data())),
        reinterpret_cast<fftw_complex *>(out.data()),
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return out;
}

std::vector<std::complex<float>> fftw(const std::vector<std::complex<float>> &in) {
    std::vector<std::complex<float>> out(in.size());
    fftwf_plan p = fftwf_plan_dft_1d(in.size(),
        reinterpret_cast<fftwf_complex *>(const_cast<std::complex<float> *>(in.data())),
        reinterpret_cast<fftwf_complex *>(out.data()),
        FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    return out;
}
#endif

#ifdef TEST_CLFFT
template<class T>
std::vector<std::complex<T>> clfft(vex::Context &ctx, const std::vector<std::complex<T>> &in) {
    typedef typename vex::cl_vector_of<T,2>::type T2;
    vex::vector<T2> input(ctx, in.size()), output(ctx, in.size());
    vex::copy(reinterpret_cast<const T2 *>(in.data()), input);
    vex::FFT<T2> fft(ctx, {in.size()}, vex::forward);
    output = fft(input);
    std::vector<std::complex<T>> out(in.size());
    vex::copy(output, reinterpret_cast<T2 *>(out.data()));
    return out;
}
#endif

#ifdef TEST_CUFFT
#include <cufft.h>
#include <cuda_runtime.h>
std::vector<std::complex<double>> cufft(const std::vector<std::complex<double>> &in) {
    size_t sz = sizeof(cufftDoubleComplex) * in.size();
    cufftHandle p; cufftPlan1d(&p, in.size(), CUFFT_Z2Z, 1);
    cufftDoubleComplex *input, *output;
    cudaMalloc(reinterpret_cast<void **>(&input), sz);
    cudaMalloc(reinterpret_cast<void **>(&output), sz);
    cudaMemcpy(input, in.data(), sz, cudaMemcpyHostToDevice);
    cufftExecZ2Z(p, input, output, CUFFT_FORWARD);
    std::vector<std::complex<double>> out(in.size());
    cudaMemcpy(out.data(), output, sz, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(p);
    cudaFree(input);
    cudaFree(output);
    return out;
}

std::vector<std::complex<float>> cufft(const std::vector<std::complex<float>> &in) {
    size_t sz = sizeof(cufftComplex) * in.size();
    cufftHandle p; cufftPlan1d(&p, in.size(), CUFFT_C2C, 1);
    cufftComplex *input, *output;
    cudaMalloc(reinterpret_cast<void **>(&input), sz);
    cudaMalloc(reinterpret_cast<void **>(&output), sz);
    cudaMemcpy(input, in.data(), sz, cudaMemcpyHostToDevice);
    cufftExecC2C(p, input, output, CUFFT_FORWARD);
    std::vector<std::complex<float>> out(in.size());
    cudaMemcpy(out.data(), output, sz, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cufftDestroy(p);
    cudaFree(input);
    cudaFree(output);
    return out;
}
#endif



int main(int argc, char **argv) {
#ifdef TEST_CLFFT
    vex::Context ctx(vex::Filter::Env && vex::Filter::Count(1));
#endif
#ifdef TEST_CUFFT
    cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
#endif

    // sizes.
    std::vector<size_t> ns;
    const size_t max_size = argc >= 3 ? boost::lexical_cast<size_t>(argv[2]) : 0;

    if(argc == 1 || std::string(argv[1]) == "all")
        for(size_t n = 1 ; n <= (max_size == 0 ? 1024 : max_size) ; n++) ns.push_back(n);
    else if(std::string(argv[1]) == "sample")
        for(size_t n = 1 ; n <= (max_size == 0 ? 1024 : max_size) ; n += rand() % 10 + 1) ns.push_back(n);
    else if(std::string(argv[1]) == "two")
        for(size_t n = 1 ; n <= (max_size == 0 ? 65536 : max_size) ; n *= 2) ns.push_back(n);
    else if(std::string(argv[1]) == "prime") {
        vex::fft::prime_generator g;
        for(size_t n = g() ; n <= (max_size == 0 ? 1024 : max_size) ; n = g()) ns.push_back(n);
    } else {
        std::cerr << "usage: " << argv[0] << " all|sample|two|prime <size> <seed>" << std::endl;
        return EXIT_FAILURE;
    }

    // data
#ifdef USE_FLOAT
    typedef float T;
#else
    typedef double T;
#endif
    typedef std::complex<T> C;
    std::vector<C> all_input(ns.back());
    std::minstd_rand gen(argc >= 4 ? boost::lexical_cast<size_t>(argv[3]) : 23);
    std::uniform_real_distribution<T> dist(-0.5, 0.5);
    for(auto &x : all_input) x = C(dist(gen), dist(gen));

    // evaluate
    std::cout << std::setprecision(10);
    std::cout << "#n\tfftw\tclfft\tcufft\tsimple\t1000" << std::endl;
    for(size_t n : ns) {
        const std::vector<C> input(all_input.begin(), all_input.begin() + n);
        const auto exact = simple_fft<float50>(input);
        std::cout << n;
#ifdef TEST_FFTW
        std::cout << '\t' << compare(exact, fftw(input));
#else
        std::cout << "\t.";
#endif
#ifdef TEST_CLFFT
        std::cout << '\t' << compare(exact, clfft(ctx, input));
#else
        std::cout << "\t.";
#endif
#ifdef TEST_CUFFT
        std::cout << '\t' << compare(exact, cufft(input));
#else
        std::cout << "\t.";
#endif
#ifdef TEST_SIMPLE
        std::cout << '\t' << compare(exact, simple_fft<T>(input));
#else
        std::cout << "\t.";
#endif
#ifdef TEST_EXACT
        std::cout << '\t' << compare(simple_fft<float1000>(input), exact);
#else
        std::cout << "\t.";
#endif
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}


