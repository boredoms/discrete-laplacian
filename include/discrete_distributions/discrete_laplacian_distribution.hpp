#ifndef _DISCRETE_LAPLACIAN_H_
#define _DISCRETE_LAPLACIAN_H_

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

  // TODO: writing out and in, this makes it conformant with the
  // RandomNumberDistribution STL trait

  // this function produces a geometrically distributed variable more quickly
  // than the standard library function, since it does not use an iterative
  // approach and continuous coin flips, but rather a single sample from the
  // uniform distribution. once p becomes large, this provides significant
  // speedups. See Knuth's AOCP Vol. 2 Section 3.4.1 Algorithm F
  template <std::uniform_random_bit_generator URNG>
  IntType fast_geometric_distribution(URNG &urng, const double p) {
    std::exponential_distribution<> ed(1);

    double u;

    do {
      u = std::generate_canonical<double, std::numeric_limits<double>::digits>(
          urng);
    } while (u == 0.0);

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
      return l._M_p == r._M_p;
    }

    friend bool operator!=(const param_type &l, const param_type &r) {
      return !(l == r);
    }
  };

  // when creating a distribution for a DP application, note that
  // $p = e^{-1/\epsillon}$
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
  result_type min() const { return std::numeric_limits<IntType>::min(); }
  result_type max() const { return std::numeric_limits<IntType>::max(); }

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

template <class CharT, class Traits, class IntegerType>
std::basic_ostream<CharT, Traits> &
operator<<(std::basic_ostream<CharT, Traits> &os,
           const DiscreteLaplacian<IntegerType> &dnd) {
  std::basic_ostream<CharT, Traits> savestate(nullptr);
  savestate.copyfmt(os);

  using OS = std::basic_ostream<CharT, Traits>;

  os.flags(OS::dec | OS::left | OS::fixed | OS::scientific);
  os << dnd.p();

  os.copyfmt(savestate);
  return os;
}

template <class CharT, class Traits, class IntegerType>
std::basic_istream<CharT, Traits> &
operator>>(std::basic_istream<CharT, Traits> &is,
           const DiscreteLaplacian<IntegerType> &dnd) {
  std::basic_istream<CharT, Traits> savestate(nullptr);
  savestate.copyfmt(is);

  using IS = std::basic_istream<CharT, Traits>;
  using param_type = DiscreteLaplacian<IntegerType>::param_type;

  is.flags(IS::dec | IS::skipws);

  double p;
  is >> p;

  if (!is.fail()) {
    dnd.param(param_type(p));
  }

  is.copyfmt(savestate);
  return is;
}

#endif // !DEBUG
