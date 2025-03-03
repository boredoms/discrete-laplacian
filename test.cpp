#include <algorithm>
#include <catch2/catch_all.hpp>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "test.hpp"

#include "canonne.hpp"
#include "discrete_gaussian_distribution.hpp"
#include "discrete_laplacian_distribution.hpp"

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
  int num_stddevs = 3;

  // set up the distribution
  DiscreteLaplacian<int> dld(p);

  auto [sample_mean, margin] =
      compute_mean_test_values(dld, num_samples, num_stddevs);

  REQUIRE_THAT(sample_mean, Catch::Matchers::WithinAbs(dld.mean(), margin));
}

// Same as the mean test but testing for the variance
TEST_CASE("Sample variance matches expected variance", "[DiscreteLaplacian]") {
  // the values of p to test with
  auto p = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);

  // how many samples are to be taken
  int num_samples = 1000000;
  int num_stddevs = 4;

  // set up the random generators
  DiscreteLaplacian<int> dld(p);

  auto [sample_variance, margin] =
      compute_variance_test_values(dld, num_samples, num_stddevs);

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
  DiscreteLaplacian<int> dld(p);

  auto [chi_2_statistic, chi_2_mean, chi_2_stddev] =
      compute_chi_test_values(dld, num_samples);

  // here we make the assumption that the chi_2 distribution is sufficiently
  // close to a normal distribution and therefore we can assume that the
  // observed statistic is usually within three times the stddev
  REQUIRE_THAT(chi_2_statistic,
               Catch::Matchers::WithinAbs(chi_2_mean, 3 * chi_2_stddev));
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
  DiscreteLaplacian<int> dld(p);

  auto buffer = generate_testing_data(dld, num_samples);

  // here we compute the test margin, we want to check all buckets within 8
  // standard distributions
  int stddev = std::ceil(std::sqrt(dld.var()));
  int margin = 8 * stddev;
  int num_buckets = 2 * margin + 1;

  auto counts = compute_counts(buffer, margin);
  int seen = std::count_if(buffer.begin(), buffer.end(),
                           [&](int x) { return x < -margin; });

  auto ks = compute_ks_statistic(dld, counts, num_samples, margin, seen);

  printf("Kolmogorov-Smirnov: %f\n", ks);

  // the numbers are from Knuth's AOCP Vol 2, and they should have the test
  // fail no more than 2% of the time
  REQUIRE((0.07089 < ks && ks < 1.5174));
}

// same as the mean test for the discrete Laplacian, but for the discrete
// Gaussian
TEST_CASE("Sample mean matches expected mean", "[DiscreteGaussian]") {
  // the values of sigma^2 to test with
  auto sigma_square = GENERATE(0.1, 0.5, 1, 2, 10);

  // how many samples are to be taken
  int num_samples = 1000000;
  int num_stddevs = 3;

  // set up the random distribution
  DiscreteGaussian<int> dnd(sigma_square);

  auto [sample_mean, margin] =
      compute_mean_test_values(dnd, num_samples, num_stddevs);

  REQUIRE_THAT(sample_mean, Catch::Matchers::WithinAbs(dnd.mean(), margin));
}

// Same as the mean test but testing for the variance
TEST_CASE("Sample variance matches expected variance", "[DiscreteGaussian]") {
  // the values of sigma^2 to test with
  auto sigma_square = GENERATE(0.1, 0.5, 1, 2, 10);

  // how many samples are to be taken
  int num_samples = 1000000;
  int num_stddevs = 4;

  DiscreteGaussian<int> dnd(sigma_square);

  auto [sample_variance, margin] =
      compute_variance_test_values(dnd, num_samples, num_stddevs);

  // as var is only an upper bound on the variance for different values of sigma
  // this is what we have to do
  REQUIRE(sample_variance < dnd.var() + margin);
}

TEST_CASE("Empirical frequency matches expected frequency",
          "[DiscreteGaussian]") {
  // the values of p to test with
  auto sigma_square = GENERATE(0.1, 0.25, 0.5, 0.75, 0.9);

  // how many samples are to be taken
  int num_samples = 1000000;

  // set up the random generators
  DiscreteGaussian<int> dnd(sigma_square);

  auto [chi_2_statistic, chi_2_mean, chi_2_stddev] =
      compute_chi_test_values(dnd, num_samples);

  // here we make the assumption that the chi_2 distribution is sufficiently
  // close to a normal distribution and therefore we can assume that the
  // observed statistic is usually within three times the stddev
  REQUIRE_THAT(chi_2_statistic,
               Catch::Matchers::WithinAbs(chi_2_mean, 3 * chi_2_stddev));
}

// here we use a chi-squared test to test whether the empirical number of
// samples is close to the expceted one for each integer within a margin around
// the mean, we test canonne's discrete laplacian
TEST_CASE("Empirical frequency matches expected frequency", "[Canonne]") {
  // the values of p to test with
  auto st = GENERATE(std::pair<int, int>(2, 5), std::pair<int, int>(3, 4),
                     std::pair<int, int>(14, 10), std::pair<int, int>(35, 10),
                     std::pair<int, int>(95, 10));
  auto [s, t] = st;

  // compute the equivalent p so we can use DiscreteLaplacian for double
  // checking
  double p = std::exp(-(static_cast<double>(s) / t));

  // how many samples are to be taken
  int num_samples = 1000000;

  // set up the random generators
  std::random_device rd;
  std::minstd_rand gen(rd());

  // generate the sample data
  auto buffer = std::vector<int>(num_samples);

  std::generate(buffer.begin(), buffer.end(),
                [&]() { return Canonne::discrete_laplace(gen, s, t); });

  DiscreteLaplacian<int> dld(p);

  // here we compute the test margin, we want to check all buckets within 8
  // standard distributions, this value is chosen arbitrarily but it covers most
  // of the outputs
  int stddev = std::ceil(std::sqrt(dld.var()));
  int margin = 8 * stddev;
  int num_buckets = 2 * margin + 1;

  auto counts = compute_counts(buffer, margin);

  // compute chi^2 statistic
  // the mean and stddev are direct properties of the chi_square distribution
  // with num_buckets degrees of freedom
  double chi_2_mean = num_buckets - 1;
  double chi_2_stddev = std::sqrt(2 * num_buckets - 2);

  auto chi_2_statistic = compute_chi_square(dld, counts, num_samples, margin);

  // here we make the assumption that the chi_2 distribution is sufficiently
  // close to a normal distribution and therefore we can assume that the
  // observed statistic is usually within three times the stddev
  REQUIRE_THAT(chi_2_statistic,
               Catch::Matchers::WithinAbs(chi_2_mean, 3 * chi_2_stddev));
}
