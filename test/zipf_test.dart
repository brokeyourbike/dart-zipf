import 'dart:math' as math;
import 'package:test/test.dart';
import 'package:zipf/zipf.dart';

/// the list of positive integers starting from 0
Iterable<int> get positiveIntegers sync* {
  var i = 0;
  while (true) yield i++;
}

void testAlpha(int n, double alpha) {
  var zipf = Zipf(n, alpha);

  // as the alpha increases, we need more samples,
  // since the frequency in the tail grows so low
  var samples = (math.pow(2.0, alpha) * 50000.0).toInt();

  var harmonic = Zipf.generalHarmonic(n, alpha);

  var buckets = <int, int>{};
  for (var i = 0; i < samples; i++) {
    var sample = zipf.sample();
    if (buckets.containsKey(sample - 1)) {
      buckets[sample - 1] += 1;
    } else {
      buckets[sample - 1] = 0;
    }
  }

  // for each bucket, see that frequency roughly matches what we expect
  // note that the ratios here are ratios _of fractions_, so 0.5 does not mean we're off by
  // 50%, it means we're off by 50% _of the frequency_. in other words, if the frequency was
  // supposed to be 0.10, and it was actually 0.05, that's a 50% difference.
  buckets.forEach((i, bucket) {
    var freq = bucket.toDouble() / samples;
    var expected = (1.0 / math.pow((i + 1).toDouble(), alpha) / harmonic);

    var offBy = (expected - freq).abs();

    // never be off by more than 10% in absolute terms
    expect(offBy, lessThan(0.1));
  });
}

void main() {
  group('Zipf', () {
    test('requires numberOfElements', () async {
      expect(() => Zipf(null, 1.0), throwsA(isA<AssertionError>()));
    });

    test('requires numberOfElements > 0', () async {
      expect(() => Zipf(-1, 1.0), throwsA(isA<AssertionError>()));
    });

    test('requires exponent', () async {
      expect(() => Zipf(1, null), throwsA(isA<AssertionError>()));
    });

    test('requires exponent > 0', () async {
      expect(() => Zipf(1, -1.0), throwsA(isA<AssertionError>()));
    });

    test('log1p', () async {
      expect(Zipf.log1p(26), 3.295836866004329);
      expect(Zipf.log1p(-1.65), isNaN);
      expect(Zipf.log1p(-1), -double.infinity);
      expect(Zipf.log1p(1.0 / 0), double.infinity);
      expect(Zipf.log1p(-0.0), -0.0);
    });

    test('expm1', () async {
      expect(Zipf.expm1(2.0), 6.38905609893065);
      expect(Zipf.expm1(-7.0), -0.9990881180344455);
      expect(Zipf.expm1(0.0), 0.0);
      expect(Zipf.expm1(1.0 / 0), double.infinity);
      expect(Zipf.expm1(-1.0 / 0), -1.0);
      expect(Zipf.expm1(0.0 / 0), isNaN);
    });

    test('valid generation', () async {
      testAlpha(10, 1.00);
      testAlpha(10, 2.00);
      testAlpha(10, 3.00);
      testAlpha(10, 1.08);
    });

    test('can get cumulativeDistribution at point', () async {
      var z = Zipf(10, 2.0);
      expect(z.cumulativeDistribution(2), 0.8065724784830178);
    });

    test('can get probability at point', () async {
      var z = Zipf(10, 2.0);
      expect(z.probability(2), 0.16131449569660355);
    });

    test('can get mean', () async {
      var z = Zipf(10, 2.0);
      expect(z.mean(), 1.8899401472010013);
    });

    test('compute the generalized harmonic number of order n of m', () async {
      expect(Zipf.generalHarmonic(10, 2.0), 1.5497677311665408);
    });
  });
}
