#ifndef _DL_TEST_H_
#define _DL_TEST_H_

#include <random>
#include <tuple>
#include <utility>
#include <vector>

// this file defines some template functions used in the tests

// helper function to generate num_samples samples from a distribution

template <typename Distribution>
std::vector<int> generate_samples(Distribution &d, std::minstd_rand &gen,
                                  int num_samples) {
  std::vector<int> buffer(num_samples);
  std::generate(buffer.begin(), buffer.end(), [&]() { return d(gen); });

  return buffer;
}

template <typename Distribution>
auto generate_testing_data(Distribution &d, int num_samples) {
std:;
  std::random_device rd;
  std::minstd_rand gen(rd());

  return generate_samples(d, gen, num_samples);
}

// helper function to calculate the sample mean
inline double calculate_sample_mean(const std::vector<int> &samples) {
  return static_cast<double>(
             std::accumulate(samples.begin(), samples.end(), 0)) /
         samples.size();
}

// helper function to compute the sample variance
inline double calculate_sample_variance(const std::vector<int> &samples) {
  return static_cast<double>(std::transform_reduce(
             samples.begin(), samples.end(), 0, std::plus<>(),
             [](double x) { return x * x; })) /
         (samples.size() - 1);
}

// helper function to compute the counts for elements within a margin from the
// mean, which is 0
inline std::vector<int> compute_counts(const std::vector<int> &buffer,
                                       int margin) {
  std::vector<int> counts(2 * margin + 1);

  for (auto n : buffer) {
    if (std::abs(n) <= margin) {
      counts[n + margin]++;
    }
  }

  return counts;
}

template <typename Distribution>
auto compute_mean_test_values(Distribution &d, int num_samples,
                              int num_stddevs) {
  auto buffer = generate_testing_data(d, num_samples);
  auto sample_mean = calculate_sample_mean(buffer);

  // here we compute the important parameters for the test, by the CLT we expect
  // the sample mean to be distributed as a normal distribution with mean 0 and
  // variance sigma/sqrt(n)
  auto stddev = std::sqrt(d.var());
  auto margin = num_stddevs * stddev / std::sqrt(num_samples);

  return std::make_pair(sample_mean, margin);
}

template <typename Distribution>
auto compute_variance_test_values(Distribution &d, int num_samples,
                                  int num_stddevs) {
  auto buffer = generate_testing_data(d, num_samples);
  auto sample_variance = calculate_sample_variance(buffer);

  // the limit stddev is the standard deviation of the random variable which
  // outputs the sample variance of our sample. this is 2 * sigma^4 / (n-1) for
  // the normal distribution that we converge towards by the central limit
  // theorem
  auto limit_stddev = std::sqrt(2 * d.var() * d.var() / (num_samples - 1));
  auto margin = 4 * limit_stddev;

  return std::make_pair(sample_variance, margin);
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

inline double chi_square_mean(int degrees_of_freedom) {
  return degrees_of_freedom - 1;
}

inline double chi_square_stddev(int degrees_of_freedom) {
  return std::sqrt(2 * static_cast<double>(degrees_of_freedom) - 2);
}

template <typename Distribution>
auto compute_chi_test_values(Distribution &d, int num_samples) {
  auto buffer = generate_testing_data(d, num_samples);

  // here we compute the test margin, we want to check all buckets within 8
  // standard distributions, this value is chosen arbitrarily but it covers most
  // of the outputs
  int stddev = std::ceil(std::sqrt(d.var()));
  int margin = 8 * stddev;
  int num_buckets = 2 * margin + 1;

  auto counts = compute_counts(buffer, margin);

  // compute chi^2 statistic
  auto chi_2 = compute_chi_square(d, counts, num_samples, margin);

  // the mean and stddev are direct properties of the chi_square distribution
  // with num_buckets degrees of freedom
  double chi_2_mean = chi_square_mean(num_buckets);
  double chi_2_stddev = chi_square_stddev(num_buckets);

  return std::make_tuple(chi_2, chi_2_mean, chi_2_stddev);
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

#endif
