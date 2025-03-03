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

// compute the Kolmogorov-Smirnov statistic, which is similar to an extended
// chi-squared test, where we check if the number of samples smaller than some
// number is what we expected
template <typename Distribution>
double compute_ks_statistic(const Distribution &dld,
                            const std::vector<int> &counts, int num_samples,
                            int margin, int seen) {
  double ks = 0;

  for (int i = -margin; i <= margin; i++) {
    auto expected = dld.cdf(i) * num_samples;
    seen += counts[i + margin];

    ks = std::max(ks, std::abs(expected - seen));
  }

  ks /= std::sqrt(num_samples);

  return ks;
}
