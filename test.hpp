#include <random>
#include <vector>

// this file defines some template functions used in the tests

// helper function to generate num_samples samples from a distribution

template <typename Distribution>
std::vector<int> generate_samples(Distribution &dld, std::minstd_rand &gen,
                                  int num_samples) {
  std::vector<int> buffer(num_samples);
  std::generate(buffer.begin(), buffer.end(), [&]() { return dld(gen); });

  return buffer;
}

// compute the chi-squared test statistic of sum((observed - expected)^2 /
// expected), which is used in the chi-square test to check if all buckets
// contain the expected number of elements
template <typename Distribution>
double compute_chi_square(const Distribution &dld,
                          const std::vector<int> &counts, int num_samples,
                          int margin) {
  double chi_2 = 0;
  for (auto i = -margin; i <= margin; i++) {
    auto expected = dld.pmf(i) * num_samples;

    chi_2 += std::pow(counts[i + margin] - expected, 2) / expected;
  }

  return chi_2;
}
