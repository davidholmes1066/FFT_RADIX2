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
 *
 *   Description: The following file contains the functions needed to
 *   compute the Radix-2 FFT. computation is done by calculating and storing
 *   repetitive values in heap memory this should make computation
 *   of the FFT more efficient, but requires extra RAM. The following
 *   functions do not use multithreading, making the algorithm suitable
 *   to be implemented on simple single core microcontrollers. Provided
 *   These at least have 16N bytes of SRAM.
 *
 *   Later additions will write data to flash memory.
 *
 *   Note: The author of this code is by no means an expert,
 *   therefore the code may not be as efficient or fast as
 *   alternatives out there. That being said,
 *   feel free to modify and improve the code to you're needs!
*/




#include "fft.h"
#include "avr_complex.h"

uint16_t calc_BitReversal(uint16_t Value)
{
    uint16_t Nr;                                                                                                 //Number of bits to be reversed
    uint16_t Temp;                                                                                               //Temporary storage
    uint16_t RValue = 0;                                                                                         //Bit reversed value of input

    Nr = (uint16_t)log2f(N)-1;                                                                                   //Computes Nr of bits to be reversed

    for(uint16_t i = 0, j = Nr; i <= Nr; i++, j--)
    {
        Temp = (0x01<<i);                                                                                        //Creates bit mask to test bit i
        RValue |= (((Value&Temp)>>i)<<j);                                                                        //Tests bit i and shifts to jth position (Bit reversal)
    }

    return RValue;                                                                                               //Returns the bit reversed Value
}



uint16_t *init_BRLookup(void)
{
    uint16_t *Lookup_Reverse = malloc(sizeof(uint16_t)*N);                                                   //Allocates memory for lookup array size(2*N)bytes

    for(uint16_t i = 0; i < N; i++)
    {
        Lookup_Reverse[i] = calc_BitReversal(i);                                                           //Calculates the bit reversal for the fft input order
    }

    return Lookup_Reverse;
}



complex float *init_WLookup(void)
{
    complex float *W = malloc(sizeof(complex float)*(N/2));                                                 //Allocate heap memory for twiddle factor lookup
    complex float TempW = 1;                                                                                     //Initial value for W0 = W^0 = 1
    complex float Wk = cexpf((-2*M_PI*I)/N);                                                                     //Value of W1 = W^1

    for(uint16_t i = 0; i < (N/2); i++)
    {
        W[i] = TempW;                                                                                            //Fills lookup array with W[k] = W^k
        TempW *= Wk;                                                                                             //Computes W^k+1
    }

    return W;                                                                                                    //Returns pointer to the memory
}



complex float *init_FFT(void)
{
    complex float* FFT_Array;
    FFT_Array = malloc(sizeof(complex float)*N);                                                            //Allocate memory for FFT
    return FFT_Array;                                                                                            //Returns pointer to the memory
}



void calc_FFT(complex float* FFT_Array, complex float* W)
{
    uint16_t PCalc = (N/2);                                                                                      //Amount of parallel butterfly computations
    complex float Temp;                                                                                          //Temporary variable for storing multiplication
    uint16_t CNr = 2;                                                                                            //Keeps track of number of calculations per step

    for(uint16_t i = 0; i < L; i++)                                                                              //Horizontal computation steps
    {
        for(uint16_t j = 0; j < PCalc; j++)                                                                      //Parallel computation steps
        {
            for(uint16_t k = 0; k < ((N/PCalc)/2); k++)                                                          //Calculation in one parallel
            {
                Temp = FFT_Array[(CNr*j)+(k+(CNr/2))] * W[k*((N/2)/(CNr/2))];                                    //Calculates multiplication in butterfly
                FFT_Array[((j*CNr)+k)+(CNr/2)] = FFT_Array[(j*CNr)+k] - Temp;                                    //Calculates and stores bottom of butterfly
                FFT_Array[(j*CNr)+k] += Temp;                                                                    //Calculates and stores top of butterfly
            }
        }

        CNr *= 2;                                                                                                //Set number of calculations per step to 2^k+1
        PCalc /= 2;                                                                                              //Set parallel computations to half
    }
}



float *init_Window(void)
{
    float *Window = malloc(sizeof(float)*(N/2));                                                            //Allocate memory for the window function
    for(uint16_t i = 0; i < (N/2); i++)
    {
        Window[i] = A0 - (A1*cosf((2*M_PI*i)/N)) + (A2*cosf((4*M_PI*i)/N)) - (A3*cosf((6*M_PI*i)/N));   //Generates 0.5Blackman-Harris window weights
    }

    return Window;
}



void apply_Window(complex float *FFT_Array, float *Window, uint16_t *Lookup_Reverse)
{
    for(uint16_t i = 0; i < (N/2); i++)
    {
        FFT_Array[Lookup_Reverse[i]] *= Window[i];                                                              //Apply window to first half of samples
    }

    for(uint16_t i = (N/2); i > 0; i--)
    {
        FFT_Array[Lookup_Reverse[N-i]] *= Window[i-1];                                                          //Apply window to second half of samples
    }
}


/*
 * The avr-gcc compiler for avr (Atmel RISC) microcontrollers does not support complex.h
 * or efficient computations with floating point numbers.
 * The functions included in avr_complex.h together with the functions below ensure fast
 * and efficient fourier transforms for the AVR microcontrollers.
 */



complexfloat *init_avr_Wlookup(void)
{
    complexfloat *W = malloc(sizeof(complexfloat)*(N/2));                                                   //Allocate heap memory for the twiddle factors
    complexfloat TempW;                                                                                          //Create temporary variable
    TempW.re = 1, TempW.im = 0;                                                                                  //Set to value W^0
    complexfloat Wk = cf_exp((-2*M_PI)/N);                                                                   //Value Wn^1

    for(uint16_t i = 0; i < (N/2); i++)
    {
        W[i].re = TempW.re;                                                                                      //Generate twiddle factors
        W[i].im = TempW.im;

        TempW = cf_multiply(TempW, Wk);                                                                    //Update temporary variable (W^(i+1))
    }

    return W;                                                                                                    //Return pointer to the complex struct containing twiddle factors
}



complexfloat  *init_avr_fft(void)
{
    complexfloat *FFT_Array = malloc(sizeof(complexfloat)*N);                                               //Allocate heap memory 4*N bytes
    return FFT_Array;                                                                                            //Returns pointer to allocated memory
}



void apply_avr_Window(complexfloat *FFT_Array, float *Window, uint16_t *Lookup_Reverse)
{
    for(uint16_t i = 0; i < (N/2); i++)
    {
        FFT_Array[Lookup_Reverse[i]] = cf_multiply_rf(FFT_Array[Lookup_Reverse[i]], Window[i]);           //Apply window to first half of samples
    }

    for(uint16_t i = (N/2); i > 0; i--)
    {
        FFT_Array[Lookup_Reverse[N-i]] = cf_multiply_rf(FFT_Array[Lookup_Reverse[N-i]], Window[i-1]);     //Apply window to second half of samples
    }
}



void calc_avr_FFT(complexfloat* FFT_Array, complexfloat* W)
{
    uint16_t PCalc = (N/2);                                                                                      //Amount of parallel butterfly computations
    complexfloat Temp;                                                                                           //Temporary variable for storing multiplication
    uint16_t CNr = 2;                                                                                            //Keeps track of number of calculations per step

    for(uint16_t i = 0; i < L; i++)                                                                              //Horizontal computation steps
    {
        for(uint16_t j = 0; j < PCalc; j++)                                                                      //Parallel computation steps
        {
            for(uint16_t k = 0; k < ((N/PCalc)/2); k++)                                                          //Calculation in one parallel
            {
                Temp = cf_multiply(FFT_Array[(CNr*j)+(k+(CNr/2))], W[k*((N/2)/(CNr/2))]);                  //Calculates multiplication in butterfly
                FFT_Array[((j*CNr)+k)+(CNr/2)] = cf_minus(FFT_Array[(j*CNr)+k],Temp);                      //Calculates and stores bottom of butterfly
                FFT_Array[(j*CNr)+k] = cf_plus(FFT_Array[(j*CNr)+k],Temp);                                 //Calculates and stores top of butterfly
            }
        }

        CNr *= 2;                                                                                                //Set number of calculations per step to 2^k+1
        PCalc /= 2;                                                                                              //Set parallel computations to half
    }
}






