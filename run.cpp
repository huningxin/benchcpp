#include "base.h"
#include "kernel-template.h"

void printResult(char *str) {
}
void printError(char *str) {
}
void printScore(char *str) {
}

int main() {
  KernelTemplate kernelTemplate;
  Base::OutputFunctions outputFunctions(printResult, printError, printScore);
  Base::benchmarks.runAll(outputFunctions, true);
}
