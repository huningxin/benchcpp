// Author: Peter Jensen

#include <stdio.h>
#include "base.h"
#include "kernel-template.h"


void printResult(char *str) {
  printf("%s", str);
}
void printError(char *str) {
  printf("%s", str);
}
void printScore(char *str) {
  printf("%s", str);
}

int main() {
  Base::OutputFunctions outputFunctions(printResult, printError, printScore);
  KernelTemplate kernelTemplate;
  Base::benchmarks.runAll(outputFunctions, true);
}
