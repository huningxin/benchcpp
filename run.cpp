// Author: Peter Jensen

#include <stdio.h>
#include "base.h"
#include "kernel-template.h"
#include "AverageFloat32x4.h"
#include "mandelbrot.h"

void printResult(char *str) {
  printf("%s\n", str);
}
void printError(char *str) {
  printf("%s\n", str);
}
void printScore(char *str) {
  printf("%s\n", str);
}

int main() {
  Base::OutputFunctions outputFunctions(printResult, printError, printScore);

  // The constructor for each of these objects will result in the benchmark being executed
  KernelTemplate    kernelTemplate;
  AverageFloat32x4  averageFloat32x4;
  Mandelbrot        mandelbrot;

  // Execute the benchmarks declared above
  Base::benchmarks.runAll(outputFunctions, true);
}
