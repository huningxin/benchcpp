#ifndef _BASE_H
#define _BASE_H

namespace Base {

typedef void (*printFunction)(char *str);

class OutputFunctions {
 public:
  OutputFunctions(printFunction printResult, printFunction printError, printFunction printScore) :
    printResult(printResult),
    printError(printError),
    printScore(printScore) {};

  printFunction printResult;
  printFunction printError;
  printFunction printScore;
};

class Kernel {
 public:
  Kernel();
};

class Benchmarks {
 public:
  static void runAll(OutputFunctions &outputFunctions, bool useAutoIterations);
  static void add(Kernel &kernel);
 private:
  
};

extern Benchmarks benchmarks;

}
#endif
