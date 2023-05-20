# dart-zipf

[![Actions Status](https://github.com/brokeyourbike/dart-zipf/workflows/build/badge.svg)](https://github.com/brokeyourbike/dart-zipf/actions?query=workflow%3Abuild)
[![codecov](https://codecov.io/gh/brokeyourbike/dart-zipf/branch/main/graph/badge.svg?token=0T3FR74Q0V)](https://codecov.io/gh/brokeyourbike/dart-zipf)

Dart implementation of a
[Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random
number generator.

This implementation is effectively a direct port of Apache Common's
[RejectionInversionZipfSampler](https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java),
written in Java. It is based on the method described by Wolfgang Hörmann and Gerhard Derflinger
in [*Rejection-inversion to generate variates from monotone discrete
distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM Transactions on Modeling
and Computer Simulation (TOMACS) 6.3 (1996)*.

Inspired by [rust-zipf](https://github.com/jonhoo/rust-zipf)

## Authors

- [Ivan Stasiuk](https://github.com/brokeyourbike) | [Twitter](https://twitter.com/brokeyourbike) | [LinkedIn](https://www.linkedin.com/in/brokeyourbike) | [stasi.uk](https://stasi.uk)

## License

[Apache-2.0e License](https://github.com/brokeyourbike/dart-zipf/blob/main/LICENSE)

