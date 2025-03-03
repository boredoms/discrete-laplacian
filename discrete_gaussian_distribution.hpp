#ifndef _DISCRETE_GAUSSIAN_H_
#define _DISCRETE_GAUSSIAN_H_

#include <cassert>
#include <cmath>
#include <numbers>
#include <random>
#include <type_traits>

#include "discrete_laplacian_distribution.hpp"

template <typename IntType = int> class DiscreteGaussian {
  static_assert(std::is_integral<IntType>::value &&
                    std::is_signed<IntType>::value,
                "template type must be a signed integral type");

public:
  using result_type = IntType;

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

  // the variance is not exactly sigma square for this distribution,
  // but sigma square provides an upper bound that is very close
  double var() const {
    if (_M_param.sigma_square() > 1.0 / 3.0) {
      return _M_param.sigma_square();
    } else {
      return 3 * std::exp(-1.0 / (2 * _M_param.sigma_square()));
    }
  }

  // this function gives an **approximate** probability mass function for the
  // discrete gaussian. it is approximate in the sense that we can not exactly
  // compute the normalization term, so we use a lower bound that is fairly
  // tight. therefore the probability mass values produced by this function
  // slightly **overestimate** the probability of the number occuring. as the
  // deviations are very small, this should not be an issue in real world
  // applications, but if the behavior of an algorithm that depends on this
  // function is not as expected, it might warrant investigation if the
  // approxmation is sufficiently accurate.
  double pmf(IntType k) {
    auto p = std::exp(-(k * k) / (2 * _M_param.sigma_square()));

    // this is an approximate normalization constant, which gives a rather tight
    // **lower bound** on the constant
    auto normalization =
        std::max(std::sqrt(2 * std::numbers::pi_v<double>) * _M_param.sigma() *
                     (1 + 2 * std::exp(-2.0 * std::numbers::pi_v<double> *
                                       std::numbers::pi_v<double> *
                                       _M_param.sigma_square())),
                 1 + 2 * std::exp(-1.0 / (2 * _M_param.sigma_square())));

    return p / normalization;
  }

  // NOTE: currently we don't implement a cumulative density function for the
  // discrete gaussian, as there is no good closed form expression

  friend bool operator==(const DiscreteGaussian &l, const DiscreteGaussian &r) {
    return l._M_param == r._M_param;
  }

  friend bool operator!=(const DiscreteGaussian &l, const DiscreteGaussian &r) {
    return !(l == r);
  }

private:
  param_type _M_param;
};

#endif
