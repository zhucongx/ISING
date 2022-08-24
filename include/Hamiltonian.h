#ifndef ISING_INCLUDE_HAMILTONIAM_H_
#define ISING_INCLUDE_HAMILTONIAM_H_
#include <random>
#include "Config.h"
class Hamiltonian {
  public:
    Hamiltonian(double j, double k);
    [[nodiscard]] double GetEnergy(const Config &config) const;
    [[nodiscard]] double GetEnergyChange(const Config &config, size_t index) const;
  private:
    double J, K;

};

#endif //ISING_INCLUDE_HAMILTONIAM_H_
