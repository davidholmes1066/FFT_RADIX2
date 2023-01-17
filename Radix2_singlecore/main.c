#include <complex.h>
#include "fft.h"

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

    for(int i = 0; i < N; i++)
    {
        FFT_Array[Lookup_Reverse[i]] = 1+1*I;                           //Test input
    }

    apply_Window(FFT_Array, Window, Lookup_Reverse);                    //Apply window function
    calc_FFT(FFT_Array, W);                                             //Calculate FFT

    for(int i = 0; i < N; i++)
    {
        printf("Index: %d\tImag: %fI\t\tReal: %f\n",i,cimagf(FFT_Array[i]),crealf(FFT_Array[i]));
    }

    free(Lookup_Reverse);                                       //free allocated array in heap
    free(Window);                                               //free allocated array in heap
    free(W);                                                    //free allocated array in heap
    free(FFT_Array);                                            //free allocated array in heap
}

