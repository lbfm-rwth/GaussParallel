[![Build Status](https://github.com/gap-packages/GaussPar/workflows/CI/badge.svg?branch=master)](https://github.com/gap-packages/GaussPar/actions?query=workflow%3ACI+branch%3Amaster)
# How to use

Note that this is not an actual gap package but is supposed to be included in the GAP package Gauss. Therefore you can't load the code by executing LoadPackage(). Instead, there are two different ways to use GaussPar.

## The Sequential Implementation

You can use this by reading `read.g` and calling Chief with IsHPC = false:
```
Read("read.g");
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := Chief(GF(q), A, numberChops, numberChops, false);;
```

`result` is a record with the components `coeffs, vectors, relations, pivotrows, pivotcols, rank`. To try it out directly type
```
result.coeffs;
```
The result's components `coeffs, vectors, relations, pivotrows, pivotcols, rank` mimic the output of the (also sequential) `EchelonMatTransformation` of the Gauss package.
For further information see `?EchelonMatTransformation`.

## The Parallel Implementation

Use the parallel implementation in a similar way. Call Chief with IsHPC = true.
Note that you need to start HPC-GAP!
```
Read("read.g");
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := Chief(GF(q), A, numberChops, numberChops, true);;
```

As in the sequential version the result is a record with the components `coeffs, vectors, relations, pivotrows, pivotcols, rank`.

## Benchmark Data

To get some benchmark data on
- how the algorithm performs in relation to the `Gauss` pkg
- wall and CPU time
- [lock contention](https://en.wikipedia.org/wiki/Lock_%28computer_science%29#Granularity)
start HPC-GAP and execute
```
Read("read.g");
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
MeasureContention(numberChops, q, A, true);
```

## Contact

Please submit bug reports, suggestions for improvements and patches via
the [issue tracker](https://github.com/lbfm-rwth/GaussPar/issues).

## License

`GaussPar` is free software you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.
