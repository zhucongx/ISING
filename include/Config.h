#ifndef ISING_ISING_INCLUDE_CONFIG_H_
#define ISING_ISING_INCLUDE_CONFIG_H_
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include "Atom.hpp"
#include "Lattice.hpp"
// using Graph = boost::adjacency_list<boost::vecS, boost::vecS>;
class Config {
  public:
    /// Constructor
    Config();
    Config(const Matrix_t &basis,
           std::vector<Lattice> lattice_vector,
           std::vector<Atom> atom_vector,
           bool update_neighbor);
    /// Getter
    [[nodiscard]] size_t GetNumAtoms() const;
    [[nodiscard]] const Matrix_t &GetBasis() const;
    [[nodiscard]] const std::vector<Lattice> &GetLatticeVector() const;
    [[nodiscard]] const std::vector<Atom> &GetAtomVector() const;
    [[nodiscard]] const std::vector<std::vector<size_t> > &GetFirstNeighborsAdjacencyList() const;
    [[nodiscard]] const std::vector<std::vector<size_t> > &GetSecondNeighborsAdjacencyList() const;
    [[nodiscard]] int GetSpinAtLatticeId(size_t lattice_id) const;
    /// Modify config
    void ChangeSpinAt(size_t lattice_id);
    /// IO
    void WriteCfg(const std::string &filename, bool neighbors_info) const;
  private:
    /// Modify config
    void ConvertRelativeToCartesian();
    void ConvertCartesianToRelative();
    void InitializeNeighborsList(size_t num_atoms);
    void UpdateNeighbors();
    /// Properties
    Matrix_t basis_{};
    std::vector<Lattice> lattice_vector_{};
    std::vector<Atom> atom_vector_{};
    // nearest neighbor lists
    std::vector<std::vector<size_t> > first_neighbors_adjacency_list_{};
    std::vector<std::vector<size_t> > second_neighbors_adjacency_list_{};
};
Vector_t GetLatticePairCenter(const Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair);
Matrix_t GetLatticePairRotationMatrix(const Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair);
std::unordered_set<size_t> GetNeighborsLatticeIdSetOfLattice(
    const Config &config, size_t lattice_id);
Config GenerateFCC(const Factor_t &factors);

#endif //ISING_ISING_INCLUDE_CONFIG_H_
