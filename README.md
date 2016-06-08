# AVX-tuned KISS FFT

this is AVX-tuned KISS FFT library.
this implementation method is different from the method in README.simd.

## Environment:
Microsoft Visual Studio 2012 Express Edition

## Build:
- debug build  
  nmake /f Makefile.vc
- release build  
  nmake /f Makefile.vc "nodebug=1"

## NOTE:
- you can use only **double** precision floating-point number for input and output.
- kissfft.hh haven't been changed to support AVX at all.

## License:
Public domain.

## AUTHOR:
Nob.Aoki  
nobaokix(at)gmail.com

## Example:
```c
#include "kiss_fft.h"
void func(kiss_fft_cpx*input, kiss_fft_cpx*output, int sz)
{
    kiss_fft_cfg pcfg = kiss_fft_alloc(sz, 0/*1:inverse*/, 0, 0);
    kiss_fft(pcfg, input, output);
    KISS_FFT_FREE(pcfg);
}
```
