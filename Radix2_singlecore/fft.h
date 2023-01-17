#ifndef RADIX2_SINGLECORE_FFT_H
#define RADIX2_SINGLECORE_FFT_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#define N 8
#define L log2f(N)

extern complex float *FFT_Array;                            //FFT array (N length) for samples & "inplace" computation
extern complex float *W;                                    //Twiddle factors (N/2 length)
extern uint16_t *Lookup_Reverse;                            //dynamic array to lookup reverse bit order

uint16_t calc_BitReversal(uint16_t Value);                  //Reverses log2(N) bits of input
void init_BRLookup(void);                                   //generates lookup array for fft input order
void init_WLookup(void);                                    //Generates lookup array for twiddle factors (Wk)
void init_FFT(void);                                        //Prepares lookup arrays and dynamically allocates memory
void calc_FFT(void);                                        //Computes radix-2 fft of N length (single core + thread)

#endif //RADIX2_SINGLECORE_FFT_H
