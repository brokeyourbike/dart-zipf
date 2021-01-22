import 'dart:math' as math;

import 'sampler.dart';

/// Random number generator that generates Zipf-distributed random numbers
/// using rejection inversion.
class Zipf extends Sampler {
  /// Creates a new [Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law)
  /// random number generator.
  ///
  /// Note that both the number of elements
  /// and the exponent must be greater than 0.
  Zipf(int numberOfElements, double exponent)
      : assert(numberOfElements != null),
        assert(exponent != null),
        assert(numberOfElements > 0),
        assert(exponent > 0.0),
        _numberOfElements = numberOfElements,
        _exponent = exponent,
        super(math.Random()) {
    _hIntegralX1 = hIntegral(1.5) - 1.0;
    _hIntegralNumberOfElements = hIntegral(_numberOfElements + 0.5);
    _s = 2.0 - hIntegralInverse(hIntegral(2.5) - h(2.0));
  }

  /// Number of elements
  final int _numberOfElements;

  /// Exponent parameter of the distribution
  final double _exponent;

  /// `hIntegral(1.5) - 1`
  double _hIntegralX1;

  /// `hIntegral(numberOfElements + 0.5)`
  double _hIntegralNumberOfElements;

  /// `hIntegralInverse(hIntegral(2.5) - h(2)`
  double _s;

  /// Rejection inversion sampling method for a discrete, bounded Zipf
  /// distribution that is based on the method described in Wolfgang Hörmann
  /// and Gerhard Derflinger.
  /// "Rejection-inversion to generate variates from monotone discrete
  /// distributions", ACM Transactions on Modeling and Computer Simulation
  /// (TOMACS) 6.3 (1996): 169-184.
  int sample() {
    // The paper describes an algorithm for exponents larger than 1
    // (Algorithm ZRI).
    // The original method uses
    //   H(x) = (v + x)^(1 - q) / (1 - q)
    // as the integral of the hat function.
    // This function is undefined for q = 1, which is the reason for
    // the limitation of the exponent.
    // If instead the integral function
    //   H(x) = ((v + x)^(1 - q) - 1) / (1 - q)
    // is used,
    // for which a meaningful limit exists for q = 1, the method works
    // for all positive exponents.
    // The following implementation uses v = 0 and generates integral
    // number in the range [1, numberOfElements].
    // This is different to the original method where v is defined to
    // be positive and numbers are taken from [0, i_max].
    // This explains why the implementation looks slightly different.

    while (true) {
      final u = _hIntegralNumberOfElements +
          nextDouble() * (_hIntegralX1 - _hIntegralNumberOfElements);
      // u is uniformly distributed in (hIntegralX1, hIntegralNumberOfElements]

      final x = hIntegralInverse(u);

      // Limit [k] to the range [1, numberOfElements] if it would be outside
      // due to numerical inaccuracies.
      var k64 = x;
      if (k64 < 1.0) {
        k64 = 1.0;
      } else if (k64 > _numberOfElements) {
        k64 = _numberOfElements.toDouble();
      }

      // float -> integer rounds towards zero
      var k = math.max(1, k64.toInt()).toInt();

      // Here, the distribution of k is given by:
      //
      //   P(k = 1) = C * (hIntegral(1.5) - hIntegralX1) = C
      //   P(k = m) = C * (hIntegral(m + 1/2) - hIntegral(m - 1/2)) for m >= 2
      //
      //   where C = 1 / (hIntegralNumberOfElements - hIntegralX1)

      if (k64 - x <= _s || u >= hIntegral(k64 + 0.5) - h(k64)) {
        // Case k = 1:
        //
        //  The right inequality is always true, because replacing k by 1 gives
        //  u >= hIntegral(1.5) - h(1) = hIntegralX1 and u is taken from
        //  (hIntegralX1, hIntegralNumberOfElements].
        //
        //  Therefore, the acceptance rate for k = 1 is P(accepted | k = 1) = 1
        //  and the probability that 1 is returned as random value is
        //  P(k = 1 and accepted) = P(accepted | k = 1) * P(k = 1) = C = C / 1^exponent
        //
        // Case k >= 2:
        //
        //  The left inequality (k - x <= s) is just a short cut
        //  to avoid the more expensive evaluation of the right inequality
        //  (u >= hIntegral(k + 0.5) - h(k)) in many cases.
        //
        //  If the left inequality is true, the right inequality is also true:
        //    Theorem 2 in the paper is valid for all positive exponents,
        //    because
        //    the requirements h'(x) = -exponent/x^(exponent + 1) < 0 and
        //    (-1/hInverse'(x))'' = (1+1/exponent) * x^(1/exponent-1) >= 0
        //    are both fulfilled.
        //    Therefore, f(x) = x - hIntegralInverse(hIntegral(x + 0.5) - h(x))
        //    is a non-decreasing function. If k - x <= s holds,
        //    k - x <= s + f(k) - f(2) is obviously also true which
        //    is equivalent to
        //    -x <= -hIntegralInverse(hIntegral(k + 0.5) - h(k)),
        //    -hIntegralInverse(u) <= -hIntegralInverse(hIntegral(k + 0.5) - h(k)),
        //    and finally u >= hIntegral(k + 0.5) - h(k).
        //
        //   Hence, the right inequality determines the acceptance rate:
        //   P(accepted | k = m) = h(m) / (hIntegrated(m+1/2) - hIntegrated(m-1/2))
        //   The probability that m is returned is given by
        //   P(k = m and accepted) = P(accepted | k = m) * P(k = m) = C * h(m) = C / m^exponent.
        //
        // In both cases the probabilities are proportional to the probability mass function
        // of the Zipf distribution.

        return k;
      }
    }
  }

  /// Computes the cumulative distribution (CDF) of the distribution at [x],
  ///  i.e. `P(X ≤ x)`.
  double cumulativeDistribution(double x) {
    if (x < 1) {
      return 0.0;
    }

    return generalHarmonic(x.toInt(), _exponent) /
        generalHarmonic(_numberOfElements, _exponent);
  }

  /// Computes the probability mass (PMF) at [k], i.e. `P(X = k)`.
  double probability(int k) {
    return (1.0 / math.pow(k, _exponent)) /
        generalHarmonic(_numberOfElements, _exponent);
  }

  /// Gets the mean of the distribution.
  double mean() {
    return generalHarmonic(_numberOfElements, _exponent - 1.0) /
        generalHarmonic(_numberOfElements, _exponent);
  }

  /// Compute the generalized harmonic number of order n of m.
  /// `(1 + 1/2^m + 1/3^m + ... + 1/n^m)`
  static double generalHarmonic(int n, double m) {
    var sum = 0.0;
    for (var i = 0; i < n; i++) {
      sum += math.pow(i + 1, -m).toDouble();
    }
    return sum;
  }

  /// `H(x)` is an integral function of `h(x)`,
  /// the derivative of `H(x)`is `h(x)`.
  double hIntegral(double x) {
    final logX = math.log(x);
    return Zipf.helper2((1 - _exponent) * logX) * logX;
  }

  /// `h(x) = 1 / x^exponent`
  double h(double x) {
    return math.exp(-_exponent * math.log(x));
  }

  /// The inverse function of `H(x)`.
  double hIntegralInverse(double x) {
    var t = x * (1.0 - _exponent);
    if (t < -1.0) {
      // Limit value to the range [-1, +inf).
      // t could be smaller than -1 in some rare cases due to numerical errors.
      t = -1.0;
    }
    return math.exp(helper1(t) * x);
  }

  /// Helper function that calculates `log(1 + x) / x`.
  /// A Taylor series expansion is used, if [x] is close to 0.
  static double helper1(final double x) {
    if (x.abs() > 1e-8) {
      return Zipf.log1p(x) / x;
    } else {
      return 1.0 - x * (0.5 - x * ((1.0 / 3.0) - 0.25 * x));
    }
  }

  /// Helper function to calculate `(exp(x) - 1) / x`.
  /// A Taylor series expansion is used, if [x] is close to 0.
  static double helper2(final double x) {
    if (x.abs() > 1e-8) {
      return Zipf.expm1(x) / x;
    } else {
      return 1.0 + x * 0.5 * (1.0 + x * (1.0 / 3.0) * (1 + 0.25 * x));
    }
  }

  /// Computes log(x + 1)
  static double log1p(double x) {
    if (x.isNaN) {
      return double.nan;
    }

    if (x.isInfinite && !x.isNegative) {
      return double.infinity;
    }

    if (x == -1.0) {
      return -double.infinity;
    }

    if (x < -1.0) {
      return double.nan;
    }

    return math.log(x + 1.0);
  }

  /// Computes exp(x) - 1
  static double expm1(double x) {
    if (x.isNaN) {
      return double.nan;
    }

    if (x.isInfinite) {
      if (x.isNegative) {
        return -1.0;
      } else {
        return double.infinity;
      }
    }

    return math.exp(x) - 1;
  }
}
