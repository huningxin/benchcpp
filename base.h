// Author: Peter Jensen

#ifndef _BASE_H
#define _BASE_H
#include <string>
#include <stdint.h>

using namespace std;

namespace Base {

class OutputFunctions;
class Configuration;
class Benchmark;
class Benchmarks;

typedef void     (*PrintFunction)(char *str);
typedef bool     (*InitFunction)(void);
typedef bool     (*CleanupFunction)(void);
typedef uint64_t (*KernelFunction)(uint64_t n);

class OutputFunctions {
 public:
  OutputFunctions(PrintFunction printResult, PrintFunction printError, PrintFunction printScore) :
    printResult(printResult),
    printError(printError),
    printScore(printScore) {};

  PrintFunction printResult;
  PrintFunction printError;
  PrintFunction printScore;
};

class Configuration {
 public:
  Configuration(string         &name,
                InitFunction    init,
                CleanupFunction cleanup,
                KernelFunction  simd,
                KernelFunction  nonSimd,
                uint64_t        iterations) :
    kernelName(name),
    kernelInit(init),
    kernelCleanup(cleanup),
    kernelSimd(simd),
    kernelNonSimd(nonSimd),
    kernelIterations(iterations) {};

  string           kernelName;
  InitFunction     kernelInit;
  CleanupFunction  kernelCleanup;
  KernelFunction   kernelSimd;
  KernelFunction   kernelNonSimd;
  uint64_t         kernelIterations;
};

class Benchmarks {
 public:
  static void runAll(OutputFunctions &outputFunctions, bool useAutoIterations);
  static void add(Benchmark *benchmark);
};

extern Benchmarks benchmarks;

class Benchmark {
 public:
  Benchmark(Configuration *config) :
    config(config),
    useAutoIterations(false),
    initOk(true),
    cleanupOk(true) {
    Base::benchmarks.add(this);
  };

  Configuration *config;
  bool           useAutoIterations;
  bool           initOk;
  bool           cleanupOk;
  uint64_t       autoIterations;
  uint64_t       actualIterations;
  uint64_t       simdTime;
  uint64_t       nonSimdTime;
};

}
#endif
