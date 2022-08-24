#include <iostream>

#include "MonteCarlo.h"
int main(int argc, char *argv[]) {
  Factor_t factor;
  if (argc != 4) {
    std::cout << "No input ." << std::endl;
    return 1;
  } else {
    for (int i = 1; i < argc; i++) {
      factor[i-1] = std::stoull(argv[i]);
    }
  }

  MonteCarlo monte_carlo(factor,
                         1e3,
                         1e4,
                         1e10);
  monte_carlo.Simulate();
  return 0;
}
