#include <utility>
#include <chrono>
#include "MonteCarlo.h"
MonteCarlo::MonteCarlo(const Factor_t &factors)
    : config_(GenerateFCC(factors)),
      log_dump_steps_(factors[0] * factors[1] * factors[2] * 4e1),
      config_dump_steps_(factors[0] * factors[1] * factors[2] * 4e3),
      maximum_number_(factors[0] * factors[1] * factors[2]* 4e7),
      early_stop_number_(factors[0] * factors[1] * factors[2]* 4e5),
      hamiltonian_(1. / 2.89, -5. / 2.89),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      index_selector_(0, config_.GetNumAtoms() - 1) {
  energy_ = hamiltonian_.GetEnergy(config_);
  std::ofstream ofs("mc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs << "factor: " << factors[0] << " " << factors[1] << " " << factors[2] << std::endl;
  ofs << "energy: " << energy_ << std::endl;
  ofs << "early_stop_number: " << early_stop_number_ << std::endl;
}
void MonteCarlo::Simulate() {
  std::ofstream ofs("mc_log.txt", std::ofstream::out | std::ofstream::app);
  ofs << "steps\tenergy_density\tenergy\tlowest_energy\tcount\n ";
  ofs.precision(8);
  auto t1 = std::chrono::high_resolution_clock::now();

  while (steps_ < maximum_number_) {

    auto lattice_id = index_selector_(generator_);
    auto dE = hamiltonian_.GetEnergyChange(config_, lattice_id);
    if (dE < 0) {
      config_.ChangeSpinAt(lattice_id);
      energy_ += dE;
    } else {
      double possibility = std::exp(-dE);
      double random_number = one_distribution_(generator_);
      if (random_number < possibility) {
        config_.ChangeSpinAt(lattice_id);
        energy_ += dE;
      }
    }

    if (energy_ < lowest_energy_ - kEpsilon) {
      count_ = 0;
    } else {
      ++count_;
    }

    Dump(ofs);
    ++steps_;

    if (count_ >= early_stop_number_) {
      break;
    }
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Simulated Annealing Monte Carlo finished in " << std::setprecision(16)
            << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds.\n";
}
void MonteCarlo::Dump(std::ofstream &ofs) {
  if (energy_ < lowest_energy_ - kEpsilon) {
    lowest_energy_ = energy_;
    // config_.WriteCfg("lowest_energy.cfg", false);
    ofs << steps_ << '\t' << energy_ / config_.GetNumAtoms() << '\t' << energy_ << '\t'
        << lowest_energy_ << '\t' << count_
        << std::endl;
  }
  if (steps_ % log_dump_steps_ == 0) {
    ofs << steps_ << '\t' << energy_ / config_.GetNumAtoms() << '\t' << energy_ << '\t'
        << lowest_energy_ << '\t' << count_
        << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) {
    config_.WriteCfg(std::to_string(steps_) + ".cfg", false);
  }
}