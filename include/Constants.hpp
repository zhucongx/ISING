#ifndef ISING_ISING_INCLUDE_CONSTANTS_HPP_
#define ISING_ISING_INCLUDE_CONSTANTS_HPP_

#include <cmath>
namespace constants {
constexpr double kLatticeConstant = 4.046;

constexpr double kFirstNearestNeighborsCutoff = 3.5;
constexpr double kSecondNearestNeighborsCutoff = 4.8;

constexpr double kNearNeighborsCutoff = 4.8;


constexpr size_t kNumThirdNearestSetSize = 60;

constexpr size_t kNumFirstNearestNeighbors = 12;
constexpr size_t kNumSecondNearestNeighbors = 6;

} // namespace constants
#endif //ISING_ISING_INCLUDE_CONSTANTS_HPP_
