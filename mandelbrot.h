// Author: Peter Jensen

#ifndef _MANDELBROT_H
#define _MANDELBROT_H

#include <stdio.h>
#include <stdint.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include "base.h"

class Mandelbrot : public Base::Benchmark {
 public:
  Mandelbrot() :
    Base::Benchmark(
      new Base::Configuration(
        string("Mandelbrot"),
        initMandelbrot,
        cleanupMandelbrot,
        simdMandelbrot,
        nonSimdMandelbrot32, // 32
        nonSimdMandelbrot64, // 64
        1000)) {}

  static uint64_t preventOptimize;

  static uint32_t mandelx132(float c_re, float c_im, uint32_t max_iterations) {
    preventOptimize++;
    float    z_re = c_re;
    float    z_im = c_im;
    uint32_t i;
    for (i = 0; i < max_iterations; ++i) {
      float z_re2 = z_re*z_re;
      float z_im2 = z_im*z_im;
      if (z_re2 + z_im2 > 4.0f) {
        break;
      }

      float new_re = z_re2 - z_im2;
      float new_im = 2.0f*z_re*z_im;
      z_re = c_re + new_re;
      z_im = c_im + new_im;
    }
    return i;
  };

  static uint32_t mandelx164(double c_re, double c_im, uint32_t max_iterations) {
    preventOptimize++;
    double    z_re = c_re;
    double    z_im = c_im;
    uint32_t i;
    for (i = 0; i < max_iterations; ++i) {
      double z_re2 = z_re*z_re;
      double z_im2 = z_im*z_im;
      if (z_re2 + z_im2 > 4.0f) {
        break;
      }

      double new_re = z_re2 - z_im2;
      double new_im = 2.0f*z_re*z_im;
      z_re = c_re + new_re;
      z_im = c_im + new_im;
    }
    return i;
  };

  static __m128i mandelx4(__m128 c_re4, __m128 c_im4, uint32_t max_iterations) {
    preventOptimize++;
    __m128   z_re4  = c_re4;
    __m128   z_im4  = c_im4;
    __m128   four4  = _mm_set_ps1(4.0f);
    __m128   two4   = _mm_set_ps1(2.0f);
    __m128i  count4 = _mm_set1_epi32(0);
    __m128i  one4   = _mm_set1_epi32(1);
    uint32_t i;

    for (i = 0; i < max_iterations; ++i) {
      __m128 z_re24 = _mm_mul_ps(z_re4, z_re4);
      __m128 z_im24 = _mm_mul_ps(z_im4, z_im4);
      __m128 mi4    = _mm_cmple_ps(_mm_add_ps(z_re24, z_im24), four4);
      if (_mm_movemask_ps(mi4) == 0) {
        break;
      }
      __m128 new_re4 = _mm_sub_ps(z_re24, z_im24);
      __m128 new_im4 = _mm_mul_ps(_mm_mul_ps(two4, z_re4), z_im4);
      z_re4 = _mm_add_ps(c_re4, new_re4);
      z_im4 = _mm_add_ps(c_im4, new_im4);
      count4 = _mm_add_epi32(count4, _mm_and_si128(_mm_castps_si128(mi4), one4));
    }
    return count4;
  };

  static bool sanityCheck() {
    uint64_t simd      = simdMandelbrot(1);
    uint64_t nonSimd32 = nonSimdMandelbrot32(1);
    uint64_t nonSimd64 = nonSimdMandelbrot64(1);
    return simd == nonSimd32 &&
           simd == nonSimd64;
  }

  static bool initMandelbrot() {
    return sanityCheck();
  }

  static bool cleanupMandelbrot() {
    return sanityCheck();
  }

  // Non SIMD versions of the kernel
  static uint64_t nonSimdMandelbrot32 (uint64_t n) {
    uint64_t result = 0;
    for (uint64_t i = 0; i < n; ++i) {
      result = mandelx132(0.01f, 0.01f, 100);
      result = mandelx132(0.01f, 0.01f, 100) | result << 8;
      result = mandelx132(0.01f, 0.01f, 100) | result << 8;
      result = mandelx132(0.01f, 0.01f, 100) | result << 8;
    }
    return result;
  }

  static uint64_t nonSimdMandelbrot64 (uint64_t n) {
    uint64_t result = 0;
    for (uint64_t i = 0; i < n; ++i) {
      result = mandelx164(0.01, 0.01, 100);
      result = mandelx164(0.01, 0.01, 100) | result << 8;
      result = mandelx164(0.01, 0.01, 100) | result << 8;
      result = mandelx164(0.01, 0.01, 100) | result << 8;
    }
    return result;
  }

  // SIMD version of the kernel
  static uint64_t simdMandelbrot (uint64_t n) {
    __m128   vec0   = _mm_set_ps1(0.01f);
    uint64_t result = 0;
    for (uint64_t i = 0; i < n; ++i) {
      __m128i r = mandelx4(vec0, vec0, 100);
      M128I_INIT(r);
      result =  (uint32_t) M128I_X(r);
      result = ((uint32_t) M128I_Y(r)) | result << 8;
      result = ((uint32_t) M128I_Z(r)) | result << 8;
      result = ((uint32_t) M128I_W(r)) | result << 8;
    }
    return result;
  }

};

uint64_t Mandelbrot::preventOptimize = 0;

#endif
