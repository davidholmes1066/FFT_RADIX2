/*
 *
 *   ██████╗  █████╗ ██████╗ ██╗██╗  ██╗     ██████╗     ███████╗███████╗████████╗
 *   ██╔══██╗██╔══██╗██╔══██╗██║╚██╗██╔╝     ╚════██╗    ██╔════╝██╔════╝╚══██╔══╝
 *   ██████╔╝███████║██║  ██║██║ ╚███╔╝█████╗ █████╔╝    █████╗  █████╗     ██║
 *   ██╔══██╗██╔══██║██║  ██║██║ ██╔██╗╚════╝██╔═══╝     ██╔══╝  ██╔══╝     ██║
 *   ██║  ██║██║  ██║██████╔╝██║██╔╝ ██╗     ███████╗    ██║     ██║        ██║
 *   ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚═╝╚═╝  ╚═╝     ╚══════╝    ╚═╝     ╚═╝        ╚═╝
 *
 *   Warning: the following FFT algorithm uses 16N bytes of heap memory!
 *
 *   Globally declare the following in the main.c file:
 *   complex float *FFT_Array;
 *   complex float *W;
 *   uint16_t *Lookup_Reverse;
 *
 *   Description: The following file contains the functions needed to
 *   compute the Radix-2 FFT. computation is done by calculating and storing
 *   repetitive values in heap memory this should make the computation
 *   of the FFT more efficient but requires more RAM. The following
 *   functions do not use multithreading, making the algorithm suitable
 *   to be implemented on simple single core microcontrollers.
 *
 *   Note: The author of this code is not an expert, therefore the code
 *   may not be as efficient or fast as alternatives out there. That
 *   being said, feel free to modify and improve the code!
*/




#include "fft.h"

uint16_t calc_BitReversal(uint16_t Value)
{
    uint16_t Nr;                                                                //Number of bits to be reversed
    uint16_t Temp;                                                              //Temporary storage
    uint16_t RValue = 0;                                                        //Bit reversed value of input

    Nr = (uint16_t)log2f(N)-1;                                                  //Computes Nr of bits to be reversed

    for(uint16_t i = 0, j = Nr; i <= Nr; i++, j--)
    {
        Temp = (0x01<<i);                                                       //Creates bit mask to test bit i
        RValue |= (((Value&Temp)>>i)<<j);                                       //Tests bit i and shifts to jth position (Bit reversal)
    }

    return RValue;                                                              //Returns the bit reversed Value
}



void init_BRLookup(void)
{
    Lookup_Reverse = malloc(sizeof(uint16_t)*N);                           //Allocates memory for lookup array size(2*N)bytes

    for(uint16_t i = 0; i < N; i++)
    {
        Lookup_Reverse[i] = calc_BitReversal(i);                          //Calculates the bit reversal for the fft input order
    }
}



void init_WLookup(void)
{
    W = malloc(sizeof(complex float)*(N/2));                               //Allocate heap memory for twiddle factor lookup
    complex float TempW = 1;                                                    //Initial value for W0 = W^0 = 1
    complex float Wk = cexpf((-2*M_PI*I)/N);                                 //Value of W1 = W^1

    for(uint16_t i = 0; i < (N/2); i++)
    {
        W[i] = TempW;                                                           //Fills lookup array with W[k] = W^k
        TempW *= Wk;                                                            //Computes W^k+1
    }
}



void init_FFT(void)
{
    init_BRLookup();                                                            //Calculate and fill reverse bit array
    init_WLookup();                                                             //Calculate and fill twiddle factor array
    FFT_Array = malloc(sizeof(complex float)*N);                           //Allocate memory for FFT
}



void calc_FFT(void)
{
    uint16_t PCalc = (N/2);                                                     //Amount of parallel butterfly computations
    complex float Temp;                                                         //Temporary variable for storing multiplication
    uint16_t CNr = 2;                                                           //Keeps track of number of calculations per step

    for(uint16_t i = 0; i < L; i++)                                             //Horizontal computation steps
    {
        for(uint16_t j = 0; j < PCalc; j++)                                     //Parallel computation steps
        {
            for(uint16_t k = 0; k < ((N/PCalc)/2); k++)                         //Calculation in one parallel
            {
                Temp = FFT_Array[(CNr*j)+(k+(CNr/2))] * W[k*((N/2)/(CNr/2))];   //Calculates multiplication in butterfly
                FFT_Array[((j*CNr)+k)+(CNr/2)] = FFT_Array[(j*CNr)+k] - Temp;   //Calculates and stores bottom of butterfly
                FFT_Array[(j*CNr)+k] += Temp;                                   //Calculates and stores top of butterfly
            }
        }

        CNr *= 2;                                                               //Set number of calculations per step to 2^k+1
        PCalc /= 2;                                                             //Set parallel computations to half
    }
}
