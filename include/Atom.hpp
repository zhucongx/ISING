#ifndef ISING_ISING_INCLUDE_ATOM_HPP_
#define ISING_ISING_INCLUDE_ATOM_HPP_
#include <fstream>
#include <map>

#include "VectorMatrix.hpp"
struct Atom {

  public:
    /// Constructor
    Atom() = default;
    Atom(size_t id, const int spin)
        : id_(id), spin_(spin) {}

    /// Getter
    [[nodiscard]] size_t GetId() const {
      return id_;
    }
    [[nodiscard]] int GetSpin() const {
      return spin_;
    }
    [[nodiscard]] std::string GetSpinString() const {
      return std::to_string(spin_);
    }
    [[nodiscard]] static constexpr double GetMass() {
      return 0;
    }
    /// Setter
    void SetSpin(const int &spin) {
      spin_ = spin;
    }
  private:
    // atom id
    size_t id_{};
    // spin type
    int spin_{};
};

#endif //ISING_ISING_INCLUDE_ATOM_HPP_
