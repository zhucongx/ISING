#ifndef ISING_INCLUDE_MONTECARLO_H_
#define ISING_INCLUDE_MONTECARLO_H_
#include <random>
#include "Hamiltonian.h"
#include "Config.h"
class MonteCarlo {
  public:
    MonteCarlo(const Factor_t &factors);
    void Simulate();
  private:
    inline void Dump(std::ofstream &ofs);
    Config config_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;
    const unsigned long long int early_stop_number_;
    // simulation statistics
    unsigned long long int steps_{0};
    unsigned long long int count_{0};
    double energy_{0.0};
    double lowest_energy_{0.0};

    const Hamiltonian hamiltonian_;
    mutable std::mt19937_64 generator_;
    mutable std::uniform_int_distribution<size_t> index_selector_;
    mutable std::uniform_real_distribution<double> one_distribution_{0.0, 1.0};
};

#endif //ISING_INCLUDE_MONTECARLO_H_
