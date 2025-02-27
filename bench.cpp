#include <catch2/catch_all.hpp>
#include <iostream>
#include <random>

#include "discrete_laplacian_distribution.hpp"

TEST_CASE("Discrete Laplacian") {
  // set up test parameters
  auto p = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);
  int num_samples = 1000000;
  std::vector<int> buffer(num_samples);

  // set up randomness
  std::random_device rd;
  std::minstd_rand gen(rd());

  DiscreteLaplacian<int> dld(p);

  std::cout << "p = " << p << "\n";

  BENCHMARK("Fill Buffer") {
    std::generate(buffer.begin(), buffer.end(), [&]() { return dld(gen); });
  };
}
