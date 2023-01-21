#include <complex.h>
#include "fft.h"
#include "avr_complex.h"

complex float *FFT_Array;                                               //FFT array (N length) for samples & "inplace" computation
complex float *W;                                                       //Twiddle factors (N/2 length)
uint16_t *Lookup_Reverse;                                               //dynamic array to lookup reverse bit order
float *Window;                                                          //Factors for window function (Blackman-Harris)

int main(void)
{
    FFT_Array = init_FFT();
    W = init_WLookup();
    Lookup_Reverse = init_BRLookup();
    Window = init_Window();

    complexfloat *FFT_avr_Array = init_avr_fft();
    complexfloat *W_avr = init_avr_Wlookup();

    for(int i = 0; i < N; i++)
    {
        FFT_Array[Lookup_Reverse[i]] = 1+1*I;                           //Test input

        FFT_avr_Array[Lookup_Reverse[i]].re = 1;                        //Test input avr complex
        FFT_avr_Array[Lookup_Reverse[i]].im = 1;
    }

    apply_Window(FFT_Array, Window, Lookup_Reverse);                    //Apply window function
    calc_FFT(FFT_Array, W);                                             //Calculate FFT

    apply_avr_Window(FFT_avr_Array,Window,Lookup_Reverse);     //AVR apply window function
    calc_avr_FFT(FFT_avr_Array,W_avr);                      //AVR apply FFT

    for(int i = 0; i < N; i++)
    {
        printf("AVRI: %f\tAVRR: %f\t\tI: %f\tR: %f\n",FFT_avr_Array[i].im, FFT_avr_Array[i].re,cimagf(FFT_Array[i]),crealf(FFT_Array[i]));
    }

    free(Lookup_Reverse);                                        //free allocated array in heap
    free(Window);                                                //free allocated array in heap
    free(W);                                                     //free allocated array in heap
    free(FFT_Array);                                             //free allocated array in heap
}

