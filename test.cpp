#include <algorithm>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numeric>
#include <random>
#include <vector>

#include "discrete_laplacian_distribution.hpp"

// helper function to generate num_samples samples from the DiscreteLaplacian
std::vector<int> generate_samples(DiscreteLaplacian<int> &dld,
                                  std::minstd_rand &gen, int num_samples) {
  std::vector<int> buffer(num_samples);
  std::generate(buffer.begin(), buffer.end(), [&]() { return dld(gen); });

  return buffer;
}

// helper function to calculate the sample mean
double calculate_sample_mean(const std::vector<int> &samples) {
  return static_cast<double>(
             std::accumulate(samples.begin(), samples.end(), 0)) /
         samples.size();
}

// helper function to compute the sample variance
double calculate_sample_variance(const std::vector<int> &samples) {
  return static_cast<double>(std::transform_reduce(
             samples.begin(), samples.end(), 0, std::plus<>(),
             [](double x) { return x * x; })) /
         (samples.size() - 1);
}

// helper function to compute the counts for elements within a margin from the
// mean, which is 0
std::vector<int> compute_counts(const std::vector<int> &buffer, int margin) {
  std::vector<int> counts(2 * margin + 1);

  for (auto n : buffer) {
    if (std::abs(n) <= margin) {
      counts[n + margin]++;
    }
  }

  return counts;
}

// compute the chi-squared test statistic of sum((observed - expected)^2 /
// expected), which is used in the chi-square test to check if all buckets
// contain the expected number of elements
double compute_chi_square(const DiscreteLaplacian<int> &dld,
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
double compute_ks_statistic(const DiscreteLaplacian<int> &dld,
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

// helper function to compute the chi_2 statistic
//

// Here we test whether the mean matches what we would expect
// Since we are summing over random variables, we can make use of the
// central limit theorem in establishing our error bounds. The variance is
// reduced by a factor of sqrt(num_samples), and we know that we should be
// within 3*std_dev 99.7% of the time
TEST_CASE("Sample mean matches expected mean", "[DiscreteLaplacian]") {
  // the values of p to test with
  auto p = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);

  // how many samples are to be taken
  int num_samples = 1000000;

  // set up the random generators and the distribution
  std::random_device rd;
  std::minstd_rand gen(rd());

  DiscreteLaplacian<int> dld(p);

  // here we compute the important parameters for the test, by the CLT we expect
  // the sample mean to be distributed as a normal distribution with mean 0 and
  // variance sigma/sqrt(n)
  auto mean = dld.mean();
  auto stddev = std::sqrt(dld.var());
  auto margin = 3 * stddev / std::sqrt(num_samples);

  // fill the buffer with random integers drawn from DiscreteLaplacian
  auto buffer = generate_samples(dld, gen, num_samples);
  auto sample_mean = calculate_sample_mean(buffer);

  REQUIRE_THAT(sample_mean, Catch::Matchers::WithinAbs(mean, margin));
}

// Same as the mean test but testing for the variance
TEST_CASE("Sample variance matches expected variance", "[DiscreteLaplacian]") {
  // the values of p to test with
  auto p = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);

  // how many samples are to be taken
  int num_samples = 1000000;

  // set up the random generators
  std::random_device rd;
  std::minstd_rand gen(rd());

  DiscreteLaplacian<int> dld(0.5);

  // the limit stddev is the standard deviation of the random variable which
  // outputs the sample variance of our sample. this is 2 * sigma^4 / (n-1) for
  // the normal distribution that we converge towards by the central limit
  // theorem
  auto limit_stddev = std::sqrt(2 * dld.var() * dld.var() / (num_samples - 1));
  auto margin = 4 * limit_stddev;

  auto buffer = generate_samples(dld, gen, num_samples);
  auto sample_variance = calculate_sample_variance(buffer);

  REQUIRE_THAT(sample_variance, Catch::Matchers::WithinAbs(dld.var(), margin));
}

// here we use a chi-squared test to test whether the empirical number of
// samples is close to the expceted one for each integer within a margin around
// the mean
TEST_CASE("Empirical frequency matches expected frequency",
          "[DiscreteLaplacian]") {
  // the values of p to test with
  auto p = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);

  // how many samples are to be taken
  int num_samples = 1000000;

  // set up the random generators
  std::random_device rd;
  std::minstd_rand gen(rd());

  DiscreteLaplacian<int> dld(0.5);

  // here we compute the test margin, we want to check all buckets within 8
  // standard distributions, this value is chosen arbitrarily but it covers most
  // of the outputs
  int stddev = std::ceil(std::sqrt(dld.var()));
  int margin = 8 * stddev;
  int num_buckets = 2 * margin + 1;

  auto buffer = generate_samples(dld, gen, num_samples);
  auto counts = compute_counts(buffer, margin);

  // compute chi^2 statistic
  // the mean and stddev are direct properties of the chi_square distribution
  // with num_buckets degrees of freedom
  double chi_2_mean = num_buckets - 1;
  double chi_2_stddev = std::sqrt(2 * num_buckets - 2);

  auto chi_2 = compute_chi_square(dld, counts, num_samples, margin);

  // here we make the assumption that the chi_2 distribution is sufficiently
  // close to a normal distribution and therefore we can assume that the
  // observed statistic is usually within three times the stddev
  REQUIRE_THAT(chi_2, Catch::Matchers::WithinAbs(chi_2_mean, 3 * chi_2_stddev));
}

// here we use the Kolmogorov-Smirnov test to test whether the empirical number
// of samples less than some integer is close to the expceted number for each
// integer within a margin around the mean
TEST_CASE("Empirical distribution matches expected distribution",
          "[DiscreteLaplacian]") {
  // the values of p to test with
  auto p = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);

  // how many samples are to be taken
  int num_samples = 1000000;

  // set up the random generators
  std::random_device rd;
  std::minstd_rand gen(rd());

  DiscreteLaplacian<int> dld(0.5);

  // here we compute the test margin, we want to check all buckets within 8
  // standard distributions
  int stddev = std::ceil(std::sqrt(dld.var()));
  int margin = 8 * stddev;
  int num_buckets = 2 * margin + 1;

  auto buffer = generate_samples(dld, gen, num_samples);

  auto counts = compute_counts(buffer, margin);
  int seen = std::count_if(buffer.begin(), buffer.end(),
                           [&](int x) { return x < -margin; });

  auto ks = compute_ks_statistic(dld, counts, num_samples, margin, seen);

  printf("Kolmogorov-Smirnov: %f", ks);

  // the numbers are from Knuth's AOCP Vol 2, and they should have the test
  // fail no more than 2% of the time
  REQUIRE((0.07089 < ks && ks < 1.5174));
}
