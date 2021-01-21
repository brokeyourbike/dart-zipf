# dart-zipf

[![Actions Status](https://github.com/flexremit/app/workflows/build/badge.svg)](https://github.com/flexremit/app/actions?query=workflow%3Abuild)
[![Codecov](https://codecov.io/github/brokeyourbike/dart-zipf/coverage.svg?branch=master)](https://codecov.io/gh/brokeyourbike/dart-zipf)
[![License: Apache-2.0](https://img.shields.io/github/license/brokeyourbike/dart-zipf)](https://github.com/brokeyourbike/dart-zipf/blob/main/LICENSE)

Dart implementation of a
[Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random
number generator.

This implementation is effectively a direct port of Apache Common's
[RejectionInversionZipfSampler](https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java),
written in Java. It is based on the method described by Wolfgang Hörmann and Gerhard Derflinger
in [*Rejection-inversion to generate variates from monotone discrete
distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM Transactions on Modeling
and Computer Simulation (TOMACS) 6.3 (1996)*.