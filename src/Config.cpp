#include "Config.h"

#include <random>
#include <chrono>
#include <utility>


Config::Config() = default;
Config::Config(const Matrix_t &basis,
               std::vector<Lattice> lattice_vector,
               std::vector<Atom> atom_vector,
               bool update_neighbor)
    : basis_(basis),
      lattice_vector_(std::move(lattice_vector)),
      atom_vector_(std::move(atom_vector)) {

  if (lattice_vector_.size() != atom_vector_.size()) {
    std::cerr << "Warning: lattice vector and atom vector size not match" << std::endl;
  }
  if (update_neighbor) {
    UpdateNeighbors();
  }
}
size_t Config::GetNumAtoms() const {
  return atom_vector_.size();
}
const Matrix_t &Config::GetBasis() const {
  return basis_;
}
const std::vector<Lattice> &Config::GetLatticeVector() const {
  return lattice_vector_;
}
const std::vector<Atom> &Config::GetAtomVector() const {
  return atom_vector_;
}
const std::vector<std::vector<size_t> > &Config::GetFirstNeighborsAdjacencyList() const {
  return first_neighbors_adjacency_list_;
}
const std::vector<std::vector<size_t> > &Config::GetSecondNeighborsAdjacencyList() const {
  return second_neighbors_adjacency_list_;
}

int Config::GetSpinAtLatticeId(size_t lattice_id) const {
  return atom_vector_[lattice_id].GetSpin();
}

void Config::ChangeSpinAt(size_t lattice_id) {
  atom_vector_.at(lattice_id).SetSpin(-GetSpinAtLatticeId(lattice_id));
}

void Config::WriteCfg(const std::string &filename, bool neighbors_info) const {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << "Number of particles = " << GetNumAtoms() << '\n';
  ofs << "A = 1.0 Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << basis_[kXDimension][kXDimension] << " A\n";
  ofs << "H0(1,2) = " << basis_[kXDimension][kYDimension] << " A\n";
  ofs << "H0(1,3) = " << basis_[kXDimension][kZDimension] << " A\n";
  ofs << "H0(2,1) = " << basis_[kYDimension][kXDimension] << " A\n";
  ofs << "H0(2,2) = " << basis_[kYDimension][kYDimension] << " A\n";
  ofs << "H0(2,3) = " << basis_[kYDimension][kZDimension] << " A\n";
  ofs << "H0(3,1) = " << basis_[kZDimension][kXDimension] << " A\n";
  ofs << "H0(3,2) = " << basis_[kZDimension][kYDimension] << " A\n";
  ofs << "H0(3,3) = " << basis_[kZDimension][kZDimension] << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom: atom_vector_) {
    size_t lattice_id = atom.GetId();
    const auto &lattice = lattice_vector_[lattice_id];
    ofs << Atom::GetMass() << '\n'
        << atom.GetSpinString() << '\n'
        << lattice.GetRelativePosition();
    if (neighbors_info) {
      ofs << " # ";
      for (auto neighbor_lattice_index: first_neighbors_adjacency_list_[lattice_id]) {
        ofs << neighbor_lattice_index << ' ';
      }
      for (auto neighbor_lattice_index: second_neighbors_adjacency_list_[lattice_id]) {
        ofs << neighbor_lattice_index << ' ';
      }
    }
    ofs << '\n';
    ofs << std::flush;
  }
}
void Config::ConvertRelativeToCartesian() {
  for (auto &lattice: lattice_vector_) {
    lattice.SetCartesianPosition(lattice.GetRelativePosition() * basis_);
  }
}
void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix(basis_);
  for (auto &lattice: lattice_vector_) {
    lattice.SetRelativePosition(lattice.GetCartesianPosition() * inverse_basis);
  }
}
void Config::InitializeNeighborsList(size_t num_atoms) {
  first_neighbors_adjacency_list_.resize(num_atoms);
  for (auto &neighbor_list: first_neighbors_adjacency_list_) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumFirstNearestNeighbors);
  }
  second_neighbors_adjacency_list_.resize(num_atoms);
  for (auto &neighbor_list: second_neighbors_adjacency_list_) {
    neighbor_list.clear();
    neighbor_list.reserve(constants::kNumSecondNearestNeighbors);
  }

}
void Config::UpdateNeighbors() {
  InitializeNeighborsList(GetNumAtoms());
  const double first_r_cutoff_square = std::pow(constants::kFirstNearestNeighborsCutoff, 2);
  const double second_r_cutoff_square = std::pow(constants::kSecondNearestNeighborsCutoff, 2);
  for (auto it1 = atom_vector_.begin(); it1 != atom_vector_.end(); ++it1) {
    for (auto it2 = atom_vector_.begin(); it2 != it1; ++it2) {
      auto first_lattice_id = it1->GetId();
      auto second_lattice_id = it2->GetId();
      Vector_t absolute_distance_vector =
          GetRelativeDistanceVectorLattice(lattice_vector_[first_lattice_id],
                                           lattice_vector_[second_lattice_id]) * basis_;
      if (std::abs(absolute_distance_vector[kXDimension])
          > constants::kNearNeighborsCutoff) { continue; }
      if (std::abs(absolute_distance_vector[kYDimension])
          > constants::kNearNeighborsCutoff) { continue; }
      if (std::abs(absolute_distance_vector[kZDimension])
          > constants::kNearNeighborsCutoff) { continue; }
      const double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square < first_r_cutoff_square) {
        first_neighbors_adjacency_list_[first_lattice_id].push_back(second_lattice_id);
        first_neighbors_adjacency_list_[second_lattice_id].push_back(first_lattice_id);
      } else if (absolute_distance_square < second_r_cutoff_square) {
        second_neighbors_adjacency_list_[first_lattice_id].push_back(second_lattice_id);
        second_neighbors_adjacency_list_[second_lattice_id].push_back(first_lattice_id);
      }
    }
  }
}

Vector_t GetLatticePairCenter(const Config &config,
                              const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  Vector_t center_position;
  for (const auto kDim: All_Dimensions) {
    double first_relative =
        config.GetLatticeVector()[lattice_id_jump_pair.first].GetRelativePosition()[kDim];
    const double second_relative =
        config.GetLatticeVector()[lattice_id_jump_pair.second].GetRelativePosition()[kDim];

    double distance = first_relative - second_relative;
    int period = static_cast<int>(distance / 0.5);
    // make sure distance is the range (0, 0.5)
    while (period != 0) {
      first_relative -= static_cast<double>(period);
      distance = first_relative - second_relative;
      period = static_cast<int>(distance / 0.5);
    }
    center_position[kDim] = 0.5 * (first_relative + second_relative);
  }
  return center_position;
}
Matrix_t GetLatticePairRotationMatrix(const Config &config,
                                      const std::pair<size_t, size_t> &lattice_id_jump_pair) {
  const auto &first_lattice = config.GetLatticeVector()[lattice_id_jump_pair.first];
  const auto &second_lattice = config.GetLatticeVector()[lattice_id_jump_pair.second];

  const Vector_t
      pair_direction = Normalize(GetRelativeDistanceVectorLattice(first_lattice, second_lattice));
  Vector_t vertical_vector{};
  for (const auto index: config.GetFirstNeighborsAdjacencyList().at(lattice_id_jump_pair.first)) {
    const Vector_t jump_vector =
        GetRelativeDistanceVectorLattice(first_lattice, config.GetLatticeVector()[index]);
    const double dot_prod = Dot(pair_direction, jump_vector);
    if (std::abs(dot_prod) < 1e-6) {
      vertical_vector = Normalize(jump_vector);
      break;
    }
  }
  // The third row is normalized since it is a cross product of two normalized vectors.
  // We use transposed matrix here because transpose of an orthogonal matrix equals its inverse
  return TransposeMatrix({pair_direction, vertical_vector,
                          Cross(pair_direction, vertical_vector)});
}

std::unordered_set<size_t> GetNeighborsLatticeIdSetOfLattice(
    const Config &config, size_t lattice_id) {
  std::unordered_set<size_t> near_neighbors_hashset;
  std::copy(config.GetFirstNeighborsAdjacencyList().at(lattice_id).begin(),
            config.GetFirstNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset,
                          near_neighbors_hashset.begin()));
  std::copy(config.GetSecondNeighborsAdjacencyList().at(lattice_id).begin(),
            config.GetSecondNeighborsAdjacencyList().at(lattice_id).end(),
            std::inserter(near_neighbors_hashset,
                          near_neighbors_hashset.begin()));

  return near_neighbors_hashset;
}
Config GenerateFCC(const Factor_t &factors) {
  static std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  static std::mt19937_64 generator(static_cast<unsigned long long int>(
                                       std::chrono::system_clock::now().time_since_epoch().count()));

  Matrix_t
      basis{{{constants::kLatticeConstant * static_cast<double>(factors[kXDimension]), 0, 0},
             {0, constants::kLatticeConstant * static_cast<double>(factors[kYDimension]), 0},
             {0, 0, constants::kLatticeConstant * static_cast<double>(factors[kZDimension])}}};
  const size_t num_atoms = 4 * factors[kXDimension] * factors[kYDimension] * factors[kZDimension];
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  std::vector<Lattice> lattice_vector;
  lattice_vector.reserve(num_atoms);
  std::vector<Atom> atom_vector;
  atom_vector.reserve(num_atoms);
  size_t count = 0;
  for (size_t k = 0; k < factors[kZDimension]; ++k) {
    for (size_t j = 0; j < factors[kYDimension]; ++j) {
      for (size_t i = 0; i < factors[kXDimension]; ++i) {
        auto x_ref = static_cast<double>(i);
        auto y_ref = static_cast<double>(j);
        auto z_ref = static_cast<double>(k);
        std::vector<Vector_t> relative_position_list = {
            {x_ref / x_length, y_ref / y_length, z_ref / z_length},
            {(x_ref + 0.5) / x_length, (y_ref + 0.5) / y_length, z_ref / z_length},
            {(x_ref + 0.5) / x_length, y_ref / y_length, (z_ref + 0.5) / z_length},
            {x_ref / x_length, (y_ref + 0.5) / y_length, (z_ref + 0.5) / z_length}
        };

        for (const auto &relative_position: relative_position_list) {
          lattice_vector.emplace_back(count, relative_position * basis, relative_position);
          atom_vector.emplace_back(count, (distribution(generator) > 0) ? 1 : -1);
          count++;
        }
      }
    }
  }
  return Config{basis, lattice_vector, atom_vector, true};
}

