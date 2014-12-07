// Author: Peter Jensen

#include <stdio.h>
#include <list>
#include <sys/types.h>
#include <sys/timeb.h>

#include "base.h"

using namespace std;

typedef list<Base::Benchmark *> BenchmarkList;

static BenchmarkList benchmarkList;

static uint64_t now() {
  struct _timeb time;
  _ftime_s(&time);
  return time.time*1000 + time.millitm;
}

static uint64_t timeKernel(Base::KernelFunction kernel, uint64_t iterations) {
  uint64_t start, stop;
  start = now();
  kernel(iterations);
  stop = now();
  return stop - start;
}

static uint64_t computeIterations(Base::Benchmark *benchmark) {
  uint64_t desiredRuntime = 1000;  // milliseconds for longest running kernel
  uint64_t testIterations = 10;    // iterations used to determine time for desiredRuntime

    // Make the slowest kernel run for at least 500ms
  uint64_t simdTime    = timeKernel(benchmark->config->kernelSimd, testIterations);
  uint64_t nonSimdTime = timeKernel(benchmark->config->kernelNonSimd, testIterations);
  uint64_t maxTime     = simdTime > nonSimdTime ? simdTime : nonSimdTime;
  while (maxTime < 500) {
    testIterations *= 2;
    simdTime = timeKernel(benchmark->config->kernelSimd, testIterations);
    nonSimdTime = timeKernel(benchmark->config->kernelNonSimd, testIterations);
    maxTime = simdTime > nonSimdTime ? simdTime : nonSimdTime;
    printf("testIterations: %llu, maxTime: %llu\n", testIterations, maxTime);
  }
  maxTime = simdTime > nonSimdTime ? simdTime : nonSimdTime;

  // Compute iteration count for 1 second run of slowest kernel
  uint64_t iterations = desiredRuntime*testIterations/maxTime;
  return iterations;
}

static bool runOne(Base::Benchmark *benchmark) {
  // Initialize the kernels and check the correctness status
  if (!benchmark->config->kernelInit()) {
    benchmark->initOk = false;
    return false;
  }

  // Determine how many iterations to use.
  if (benchmark->useAutoIterations) {
    benchmark->autoIterations = computeIterations(benchmark);
    benchmark->actualIterations = benchmark->autoIterations;
  }
  else {
    benchmark->actualIterations = benchmark->config->kernelIterations;
  }

  // Run the SIMD kernel
  benchmark->simdTime = timeKernel(benchmark->config->kernelSimd, benchmark->actualIterations);

  // Run the non-SIMD kernel
  benchmark->nonSimdTime = timeKernel(benchmark->config->kernelNonSimd, benchmark->actualIterations);

  // Do the final sanity check
  if (!benchmark->config->kernelCleanup()) {
    benchmark->cleanupOk = false;
    return false;
  }
  return true;
}

static void report(Base::Benchmark *benchmark, Base::OutputFunctions &outputFunctions) {
  char buf[200];
  if (!benchmark->initOk) {
    sprintf_s(buf, "%s: %s", benchmark->config->kernelName.c_str(), "FAILED INIT");
    outputFunctions.printError(buf);
    return;
  }
  if (!benchmark->cleanupOk) {
    sprintf_s(buf, "%s: %s", benchmark->config->kernelName.c_str(), "FAILED CLEANUP");
    outputFunctions.printError(buf);
    return;
  }
  double ratio = (double)benchmark->nonSimdTime / (double)benchmark->simdTime;
  sprintf_s(buf, "%-23s : Iterations(%10llu), SIMD(%8llums), Non-SIMD(%8llums), Speedup(%.3f)",
            benchmark->config->kernelName.c_str(),
            benchmark->actualIterations,
            benchmark->simdTime, 
            benchmark->nonSimdTime, ratio);
  outputFunctions.printResult(buf);
}

void Base::Benchmarks::runAll(Base::OutputFunctions &outputFunctions, bool useAutoIterations) {
//  printf("runAll\n");
  for (BenchmarkList::iterator it = benchmarkList.begin(); it != benchmarkList.end(); it++) {
    Base::Benchmark *benchmark = *it;
//    printf("Running: %s\n", benchmark->config->kernelName.c_str());
    benchmark->useAutoIterations = useAutoIterations;
    runOne(benchmark);
    report(benchmark, outputFunctions);
  }
}

void Base::Benchmarks::add(Base::Benchmark *benchmark) {
//  printf("adding: %s\n", benchmark->config->kernelName.c_str());
  benchmarkList.push_back(benchmark);
}
