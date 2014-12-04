#include <stdio.h>
#include "base.h"

Base::Kernel::Kernel() {
  Base::benchmarks.add(*this);
};

void Base::Benchmarks::runAll(Base::OutputFunctions &outputFunctions, bool useAutoIterations) {
  printf("runAll\n");
}

void Base::Benchmarks::add(Base::Kernel &kernel) {
  printf("add\n");
}
