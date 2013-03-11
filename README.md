`multi_precision_fft.hpp` implements a simple FFT for accuracy tests, based on the algorithms used in FFTw's [benchFFT](http://www.fftw.org/benchfft/).
Compatible with native types and [boost::multiprecision](http://www.boost.org/libs/multiprecision/), easier to adjust to any precision than the original code.

`bench.cpp` compares the precision of [FFTw](http://www.fftw.org/), [NVidia cuFFT](https://developer.nvidia.com/cufft) and [VexCL's FFT](https://ddemidov.github.com/vexcl/), but can be trivially adapted.

----------------

This software is in the public domain, furnished "as is", without technical
support, and with no warranty, express or implied, as to its usefulness for
any purpose.
