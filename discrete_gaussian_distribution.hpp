#include <cassert>
#include <cmath>
#include <random>
#include <type_traits>

#include "discrete_laplacian_distribution.hpp"

template <typename IntType = int> class DiscreteGaussian {
  static_assert(std::is_integral<IntType>::value &&
                    std::is_signed<IntType>::value,
                "template type must be a signed integral type");

public:
  using result_type = IntType;

  // TODO: implement
  class param_type {
    double _M_sigma_square;
    double _M_sigma;

  public:
    explicit param_type(double sigma_square)
        : _M_sigma_square(sigma_square), _M_sigma(std::sqrt(sigma_square)) {
      assert(_M_sigma_square > 0 && "the variance must be positive");
    }

    double sigma_square() const { return _M_sigma_square; }
    double sigma() const { return _M_sigma; }

    friend bool operator==(const param_type &l, const param_type &r) {
      return l._M_sigma_square == r._M_sigma_square;
    }

    friend bool operator!=(const param_type &l, const param_type &r) {
      return !(l == r);
    }
  };

  explicit DiscreteGaussian(double mu, double sigma) : _M_param(mu, sigma) {}
  DiscreteGaussian(const param_type &param) : _M_param(param) {}

  void reset() {}

  template <std::uniform_random_bit_generator URNG>
  result_type operator()(URNG &urng) {
    return (*this)(urng, _M_param);
  }

  // TODO: implement
  template <std::uniform_random_bit_generator URNG>
  result_type operator()(URNG &urng, const param_type &param) {
    auto t = std::floor(param.sigma()) + 1;

    while (1) {
      // sample discrete laplacian
      auto y = DiscreteLaplacian<>(t)(urng);
      // sample bernoulli
      auto p = std::exp(-std::pow(std::abs(y) - param.sigma_square() / t, 2) /
                        (2 * param.sigma_square()));
      auto c = std::bernoulli_distribution(p)(urng);

      if (c == 1) {
        return y;
      }
    }
  }

  double sigma_square() const { return _M_param.sigma_square(); }
  param_type param() const { return _M_param; }

  // these are defined by the integer type we use
  result_type min() const { std::numeric_limits<IntType>::min(); }
  result_type max() const { std::numeric_limits<IntType>::max(); }

  // the discrete gaussian is symmetric around 0
  result_type mean() const { return 0; }

private:
  param_type _M_param;
};
