#include "Hamiltonian.h"
Hamiltonian::Hamiltonian(double j, double k) : J(j), K(k) {}
double Hamiltonian::GetEnergy(const Config &config) const {
  double energy = 0.0;
  const auto &first_neighbors_adjacency_lists = config.GetFirstNeighborsAdjacencyList();
  const auto &second_neighbors_adjacency_lists = config.GetSecondNeighborsAdjacencyList();
  for (size_t i = 0; i < config.GetNumAtoms(); ++i) {
    int spin1 = config.GetAtomVector()[i].GetSpin();
    for (auto j: first_neighbors_adjacency_lists[i]) {
      int spin2 = config.GetAtomVector()[j].GetSpin();
      energy += J * spin1 * spin2;
    }
  }
  for (size_t i = 0; i < config.GetNumAtoms(); ++i) {
    int spin1 = config.GetAtomVector()[i].GetSpin();
    for (auto j: second_neighbors_adjacency_lists[i]) {
      int spin2 = config.GetAtomVector()[j].GetSpin();
      energy += K * spin1 * spin2;
    }
  }
  return energy / 2;
}
double Hamiltonian::GetEnergyChange(const Config &config, size_t index) const {
  double energy_change = 0.0;
  int old_spin = config.GetAtomVector()[index].GetSpin();
  int new_spin = -old_spin;

  for (auto j: config.GetFirstNeighborsAdjacencyList()[index]) {
    int spin2 = config.GetAtomVector()[j].GetSpin();
    energy_change += J * new_spin * spin2;
    energy_change -= J * old_spin * spin2;
  }
  for (auto j: config.GetSecondNeighborsAdjacencyList()[index]) {
    int spin2 = config.GetAtomVector()[j].GetSpin();
    energy_change += J * new_spin * spin2;
    energy_change -= J * old_spin * spin2;
  }

  return energy_change;
}
