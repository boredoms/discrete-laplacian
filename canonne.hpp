#include <cassert>
#include <random>

// NOTE: this code is provided for completeness, i do not recommend using it as
// it performs worse than the other algorithms for no apparent payoff.

// Algorithm 1 from the Canonne paper to sample from Bernoulli(exp(-gamma))
template <std::uniform_random_bit_generator URNG>
bool bernoulli_exp(URNG &urng, double gamma) {
  assert(gamma >= 0.0 && "Gamma must be non-negative");

  if (gamma <= 1.0) {
    int k = 1;

    while (true) {
      std::bernoulli_distribution bd(gamma / k);
      auto a = bd(urng);

      if (a == 0) {
        break;
      } else {
        k++;
      }
    }

    return k % 2;
  } else {
    std::bernoulli_distribution bd(std::exp(-1));

    for (int k = 1; k <= static_cast<int>(floor(gamma)); k++) {
      auto b = bd(urng);

      if (b == 0) {
        return 0;
      }
    }
    return bernoulli_exp(urng, gamma - std::floor(gamma));
  }
}

// an algorithm to sample from DiscreteLaplacian(exp^(-s/t)). This is
// algorithm 2 in the Canonne paper.
template <std::uniform_random_bit_generator URNG>
int discrete_laplace(URNG &urng, int s, int t) {
  std::uniform_int_distribution<int> ud(0, t - 1);

  while (true) {
    auto u = ud(urng);
    auto d = bernoulli_exp(urng, static_cast<double>(u) / t);

    if (d == 0) {
      continue;
    }

    int v = 0;
    std::bernoulli_distribution bd(std::exp(-1));

    while (true) {
      auto a = bd(urng);
      if (a == 0) {
        break;
      } else {
        v++;
      }
    }

    double x = u + t * v;
    auto y = std::floor(x / s);

    auto b = std::bernoulli_distribution()(urng);

    if (b == 1 && y == 0) {
      continue;
    } else {
      return (1 - 2 * b) * static_cast<int>(y);
    }
  }
}
