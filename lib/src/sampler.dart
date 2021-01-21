import 'dart:math' as math;

/// Base class for a sampler.
class Sampler {
  /// Generator of uniformly distributed random numbers.
  Sampler(math.Random random) : _random = random;

  final math.Random _random;

  /// Generates a non-negative random floating point value uniformly
  /// distributed in the range from 0.0, inclusive, to 1.0, exclusive.
  double nextDouble() {
    return _random.nextDouble();
  }

  /// Generates a non-negative random integer uniformly distributed
  /// in the range from 0, inclusive, to max, exclusive.
  int nextInt(int max) {
    return _random.nextInt(max);
  }
}
