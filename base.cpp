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
//    printf("testIterations: %llu, maxTime: %llu\n", testIterations, maxTime);
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

static void printHeaders(Base::PrintFunction printFunction) {
  char buf[200];
  sprintf_s(buf, "%-20s : %12s %12s %12s %12s %10s %10s",
            "Name", "Iterations", "Scalar32(ns)", "Scalar64(ns)", "SIMD32(ns)", "Ratio32", "Ratio64");
  printFunction(buf);
}

static void printColumns(Base::PrintFunction printFunction, const char *name, uint64_t iterations, uint64_t scalar32, uint64_t scalar64, uint64_t simd32, double ratio32, double ratio64) {
  char buf[200];
  sprintf_s(buf, "%-20s : %12llu %12llu %12llu %12llu %10.2f %10.2f",
            name, iterations, scalar32, scalar64, simd32, ratio32, ratio64);
  printFunction(buf);
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
  double ratio32 = (double)benchmark->nonSimdTime / (double)benchmark->simdTime;
  double ratio64 = (double)benchmark->nonSimdTime / (double)benchmark->simdTime;
  printColumns(
    outputFunctions.printResult,
    benchmark->config->kernelName.c_str(),
    benchmark->actualIterations,
    benchmark->nonSimdTime*1000*1000/benchmark->actualIterations,
    benchmark->nonSimdTime*1000*1000/benchmark->actualIterations,
    benchmark->simdTime*1000*1000/benchmark->actualIterations,
    ratio32,
    ratio64);
}

void Base::Benchmarks::runAll(Base::OutputFunctions &outputFunctions, bool useAutoIterations) {
//  printf("runAll\n");
  printHeaders(outputFunctions.printResult);
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
