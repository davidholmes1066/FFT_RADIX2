#ifndef RADIX2_SINGLECORE_FFT_H
#define RADIX2_SINGLECORE_FFT_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "avr_complex.h"

#define N 8                                                                                 //Defines the FFT bins
#define L log2f(N)

#define A0 0.35875                                                                          //Factors Blackman-Harris window
#define A1 0.48829
#define A2 0.14128
#define A3 0.01168



uint16_t calc_BitReversal(uint16_t Value);                                                  //Reverses log2(N) bits of input
uint16_t *init_BRLookup(void);                                                              //generates lookup array for fft input order (returns pointer to memory)
complex float *init_WLookup(void);                                                          //Generates lookup array for twiddle factors (Wk) (returns pointer to memory)
complex float *init_FFT(void);                                                              //Allocates FFT memory for "in place" computation (returns pointer to memory)
void calc_FFT(complex float *FFT_Array, complex float *W);                                  //Computes radix-2 fft of N length (single core + thread) (takes pointer to FFT_Array and W)
float *init_Window(void);                                                                   //Generates lookup array with window factors
void apply_Window(complex float *FFT_Array, float *Window, uint16_t *Lookup_Reverse);       //Applies window function to the input samples

complexfloat *init_avr_Wlookup(void);                                                       //Avr version of init_BRLookup
complexfloat  *init_avr_fft(void);                                                          //Avr version of init_fft
void apply_avr_Window(complexfloat *FFT_Array, float *Window, uint16_t *Lookup_Reverse);    //Avr version of apply_Window
void calc_avr_FFT(complexfloat* FFT_Array, complexfloat* W);                                //Avr version of calc_FFT


#endif //RADIX2_SINGLECORE_FFT_H
