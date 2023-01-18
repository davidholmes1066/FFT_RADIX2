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

    complexfloat E;
    complex float T;

    E = cf_exp((float)(2*M_PI/N));
    T = cexpf(2*I*M_PI/N);

    printf("E: I: %f\tR: %f\n",E.im, E.re);
    printf("T: I: %f\tR: %f\n",cimagf(T), crealf(T));

    for(int i = 0; i < N; i++)
    {
        FFT_Array[Lookup_Reverse[i]] = 1+1*I;                           //Test input
    }

    apply_Window(FFT_Array, Window, Lookup_Reverse);                    //Apply window function
    calc_FFT(FFT_Array, W);                                             //Calculate FFT

    free(Lookup_Reverse);                                        //free allocated array in heap
    free(Window);                                                //free allocated array in heap
    free(W);                                                     //free allocated array in heap
    free(FFT_Array);                                             //free allocated array in heap
}

