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
complexfloat cf_exp(float phi);                                                     //Computes complex exponent of phi

#endif //RADIX2_SINGLECORE_AVR_COMPLEX_H
