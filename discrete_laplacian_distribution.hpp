#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <random>
#include <type_traits>

template <typename IntType = int> class DiscreteLaplacian {
  static_assert(std::is_integral<IntType>::value &&
                    std::is_signed<IntType>::value,
                "template type must be a signed integral type");

  // TODO: profile and see if own geometric distribution based on truncation of
  // exponential distribution is faster
  //
  // TODO: writing out and in

  template <std::uniform_random_bit_generator URNG>
  IntType fast_geometric_distribution(URNG &urng, const double p) {
    double u = 0.0;

    while (u == 0.0) {
      u = std::generate_canonical<double, 53>(urng);
    }

    return static_cast<IntType>(std::ceil(std::log(u) / std::log(1 - p)));
  }

public:
  using result_type = IntType;
  using result_float_type = double;

  class param_type {
    double _M_p;

  public:
    using distribution_type = DiscreteLaplacian;

    explicit param_type(double p) : _M_p(p) {
      assert((_M_p > 0.0) && (_M_p < 1.0) && "p must be between 0 and 1");
    }

    double p() const { return _M_p; }

    friend bool operator==(const param_type &l, const param_type &r) {
      l._M_p == r._M_p;
    }

    friend bool operator!=(const param_type &l, const param_type &r) {
      return !(l == r);
    }
  };

  explicit DiscreteLaplacian(double p) : _M_param(p) {}
  DiscreteLaplacian(const param_type &param) : _M_param(param) {}

  void reset() {}

  template <std::uniform_random_bit_generator URNG>
  result_type operator()(URNG &urng) {
    return (*this)(urng, _M_param);
  }

  template <std::uniform_random_bit_generator URNG>
  result_type operator()(URNG &urng, const param_type &param) {
    return fast_geometric_distribution(urng, 1 - param.p()) -
           fast_geometric_distribution(urng, 1 - param.p());
  }

  double p() const { return _M_param.p(); }
  param_type param() const { return _M_param; }

  // these are defined by the integer type we use
  result_type min() const { std::numeric_limits<IntType>::min(); }
  result_type max() const { std::numeric_limits<IntType>::max(); }

  // the mean (first moment) of the distribution is 0, c.f. Proposition 2.2. in
  // Inusah and Kozubowski
  result_type mean() const { return 0; }

  // the second moment of the distribution is 2p/(1 − p)^2, c.f.
  // Proposition 2.2. in Inusah and Kozubowski
  result_float_type var() const {
    auto p = _M_param.p();
    auto q = 1 - p;

    return 2 * p / (q * q);
  }

  // calculate the probability mass function for the integer k
  result_float_type pmf(IntType k) const {
    auto p = _M_param.p();

    return (1 - p) / (1 + p) * std::pow(p, std::abs(k));
  }

  // calculate the cumulative probability density function P(DL(p) <= k)
  result_float_type cdf(IntType k) const {
    auto p = _M_param.p();

    if (k <= 0) {
      return std::pow(p, std::ceil(-k)) / (1 + p);
    } else {
      return 1 - std::pow(p, std::ceil(k) + 1) / (1 + p);
    }
  }

  friend bool operator==(const DiscreteLaplacian &l,
                         const DiscreteLaplacian &r) {
    return l._M_param == r._M_param;
  }

  friend bool operator!=(const DiscreteLaplacian &l,
                         const DiscreteLaplacian &r) {
    return !(l == r);
  }

private:
  param_type _M_param;
};
