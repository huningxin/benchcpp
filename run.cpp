// Author: Peter Jensen

#include <stdio.h>
#include "base.h"
#include "kernel-template.h"
#include "AverageFloat32x4.h"

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
  KernelTemplate   kernelTemplate;
  AverageFloat32x4 averageFloat32x4;
  Base::benchmarks.runAll(outputFunctions, true);
}
