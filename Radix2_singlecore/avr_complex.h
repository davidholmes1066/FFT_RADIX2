//
// Created by Dell Latitude on 18/01/2023.
//

#ifndef RADIX2_SINGLECORE_AVR_COMPLEX_H
#define RADIX2_SINGLECORE_AVR_COMPLEX_H

#include <math.h>

typedef struct
{
    float im;
    float re;
} complexfloat;

complexfloat cf_multiply(complexfloat A, complexfloat B);                           //Multiplication of two complex numbers
complexfloat cf_multiply_rf(complexfloat A, float B);                               //Multiplication float by real number
complexfloat cf_exp(float phi);                                                     //Computes complex exponent of phi (cf_exp(phi) == cexpf(phi*I))
complexfloat cf_plus(complexfloat A, complexfloat B);                               // complex A + complex B
complexfloat cf_minus(complexfloat A, complexfloat B);                              //complex A - complex B

#endif //RADIX2_SINGLECORE_AVR_COMPLEX_H
