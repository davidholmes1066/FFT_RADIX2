cmake_minimum_required(VERSION 3.24)
project(Radix2_singlecore C)

set(CMAKE_C_STANDARD 99)

add_executable(Radix2_singlecore main.c fft.c fft.h avr_complex.c avr_complex.h)
