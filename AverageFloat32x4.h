// Author: Peter Jensen

#ifndef _AVERAGEFLOAT32X4_H
#define _AVERAGEFLOAT32X4_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "base.h"

class AverageFloat32x4 : public Base::Benchmark {
 public:
  AverageFloat32x4() :
    Base::Benchmark(
      new Base::Configuration(
        string("AverageFloat32x4"),
        initArray,
        cleanup,
        simdAverage,
        average32,
        average64,
        1000)) {}

  static uint64_t preventOptimize;

  static const uint32_t length = 10000;
  static float a[length];

  static bool sanityCheck() {
    float simdVal      = simdAverageKernel();
    float nonSimd32Val = nonSimdAverageKernel32();
    float nonSimd64Val = (float) nonSimdAverageKernel64();
    return fabs(simdVal - nonSimd32Val) < 0.0001 &&
           fabs(simdVal - nonSimd64Val) < 0.0001;
  }

  static bool initArray() {
    for (uint32_t i = 0; i < length; ++i) {
      a[i] = 0.1f;
    }
    // Check that the two kernel functions yields the same result.
    return sanityCheck();
  }

  static bool cleanup() {
    for (uint32_t i = 0; i < length; ++i) {
      a[i] = 0.1f;
    }
    return sanityCheck();
  };

  static float simdAverageKernel() {
    preventOptimize++;
    __m128 sumx4 = _mm_set_ps1(0.0);
    for (uint32_t j = 0, l = length; j < l; j = j + 4) {
      sumx4 = _mm_add_ps(sumx4, _mm_loadu_ps(&(a[j])));
    }
    Base::Lanes<__m128, float> lanes(sumx4);
    return (lanes.x() + lanes.y() + lanes.z() + lanes.w())/length;
//    M128_INIT(sumx4);
//    return (M128_X(sumx4) + M128_Y(sumx4) + M128_Z(sumx4) + M128_W(sumx4))/length;
  }

  static float nonSimdAverageKernel32() {
    preventOptimize++;
    float sum = 0.0;
    for (uint32_t j = 0, l = length; j < l; ++j) {
      sum += a[j];
    }
    return sum/length;
  }

  static double nonSimdAverageKernel64() {
    preventOptimize++;
    double sum = 0.0;
    for (uint32_t j = 0, l = length; j < l; ++j) {
      sum += (double)a[j];
    }
    return sum/length;
  }

  static uint64_t simdAverage(uint64_t n) {
    float val;
    for (uint64_t i = 0; i < n; ++i) {
      val = simdAverageKernel();
    }
    return (uint64_t)val;
  };

  static uint64_t average32(uint64_t n) {
    float val;
    for (uint64_t i = 0; i < n; ++i) {
      val = nonSimdAverageKernel32();
    }
    return (uint64_t)val;
  };

  static uint64_t average64(uint64_t n) {
    double val;
    for (uint64_t i = 0; i < n; ++i) {
      val = nonSimdAverageKernel64();
    }
    return (uint64_t)val;
  };

};

uint64_t AverageFloat32x4::preventOptimize = 0;
float AverageFloat32x4::a[AverageFloat32x4::length];

#endif
