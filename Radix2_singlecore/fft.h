#ifndef RADIX2_SINGLECORE_FFT_H
#define RADIX2_SINGLECORE_FFT_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#define N 8                                                     //Defines the FFT bins
#define L log2f(N)


uint16_t calc_BitReversal(uint16_t Value);                      //Reverses log2(N) bits of input
uint16_t* init_BRLookup(void);                                  //generates lookup array for fft input order (returns pointer to memory)
complex float* init_WLookup(void);                              //Generates lookup array for twiddle factors (Wk) (returns pointer to memory)
complex float* init_FFT(void);                                  //Allocates FFT memory for "in place" computation (returns pointer to memory)
void calc_FFT(complex float* FFT_Array, complex float* W);      //Computes radix-2 fft of N length (single core + thread) (takes pointer to FFT_Array and W)

#endif //RADIX2_SINGLECORE_FFT_H
