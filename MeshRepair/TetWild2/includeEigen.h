#pragma once

// We require at least SSE 4.1, also SSSE3
#define EIGEN_VECTORIZE_SSE3
#define EIGEN_VECTORIZE_SSSE3
#define EIGEN_VECTORIZE_SSE4_1

#ifdef __AVX2__
// All Intel and AMD processors which support AVX2 also support FMA3
// The only exception is VIA Nano QuadCore: https://en.wikipedia.org/wiki/List_of_VIA_Nano_microprocessors#Nano_C
// If that's your target platform, simply remove the following line:
#define EIGEN_VECTORIZE_FMA
#endif

// Quite a few places in this project use Eigen's dense matrices (Eigen's vector type is a typedef), and related algorithms
#include <Eigen/Dense>