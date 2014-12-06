// Author: Peter Jensen

#ifndef _KERNEL_TEMPLATE_H
#define _KERNEL_TEMPLATE_H

#include <stdio.h>
#include "base.h"

class KernelTemplate : public Base::Benchmark {
 public:
  KernelTemplate() :
    Base::Benchmark(
      new Base::Configuration(
        string("kernel template"),
        init,
        cleanup,
        simd,
        nonSimd,
        1000)) {}

  static bool init() {
    return true;
  };

  static bool cleanup() {
    return true;
  };

  static uint64_t simd(uint64_t n) {
    uint64_t s = 0;
    for (uint64_t i = 0; i < n; ++i) {
      s += i;
    }
    return s;
  };

  static uint64_t nonSimd(uint64_t n) {
    uint64_t s = 0;
    for (uint64_t i = 0; i < n; ++i) {
      s += i;
    }
    return s;
  };

};
#endif
