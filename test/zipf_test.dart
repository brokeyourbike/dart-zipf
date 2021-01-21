import 'package:test/test.dart';
import 'package:zipf/zipf.dart';

void main() {
  group('Zipf', () {
    test('requires numberOfElements', () async {
      expect(() => Zipf(null, 1.0), throwsA(isA<AssertionError>()));
    });

    test('requires exponent', () async {
      expect(() => Zipf(1, null), throwsA(isA<AssertionError>()));
    });
  });
}
