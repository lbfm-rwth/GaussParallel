# Introduction

This is an implementation of a parallel Gaussian algorithm as described in
[A parallel algorithm for Gaussian elemination over finite fields](
https://arxiv.org/abs/1806.04211),
which divides the input matrix into blocks and the Gaussian algorithm into
parts which can be run in parallel.

It can be used in GAP, running on only one core, and in HPC-GAP.

# Documentation

Call makedoc.g to create the documentation.

# Installation

Download the `GaussPar` package from [here](
https://github.com/lbfm-rwth/GaussPar/archive/master.zip) and extract it into
the `pkg/` folder of your GAP or HPC-GAP installation.
To load the package open a GAP session and type:
```
LoadPackage("GaussPar");
```

# Usage

Until now there are three different functions that you can use. The function
`EchelonMatBlockwise` and `EchelonMatTransformationBlockwise` resemble
`EchelonMat` and `EchelonMatTransformation` from the Gauss package. They need a
matrix as input and return a record with different amounts of information.
The third function is called `DoEchelonMatTransformationBlockwise`. As input
you can specify the matrix of course, the field of the matrix, the numbers of
blocks that are used and whether to use the sequential or the parallel version.

## EchelonMatBlockwise

The most basic usage of our algorithm would be the following:
```
n := 4000;; numberBlocks := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := EchelonMatBlockwise(A);;
```

`result` is a record with the components `vectors` and `heads`, whereas
`vectors` is a list of the nonzero rows of the reduced echelon form of the
matrix and `heads` is a list that of the pivot elements. If the i-th entry of
`heads` is j, then we now that there is a pivot element at the position [j,i].

## EchelonMatTransformationBlockwise

To receive more information you should to this:
```
n := 4000;; numberBlocks := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := EchelonMatTransformationBlockwise(A);;
```

Now the record `result` does not contain only `vectors` and `heads`, but also
`coeffs` and `relations`, whereas `coeffs` is the transformation matrix that
the algorithm used to obtain the reduced row echelon form of the matrix and
`relations` is the kernel of the matrix.

## DoEchelonMatTransformationBlockwise

There is another function that you can use when you want to specify a few
things yourself. In a record you can specify `galoisField`,
`numberBlocksHeight`, `numberBlocksWidth`, `withTrafo` and `verify`. Whereas
`withTrafo := true` makes the function return the transformation matrix and
`verify := true` verifies the result via the transformation matrix.

```
n := 4000;; numberBlocks := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
# If you want to specify the field
result := DoEchelonMatTransformationBlockwise(A, rec( galoisField := GF(q) ));;

# numberBlocksHeight and numberBlocksWidth let you determine what number
# of blocks there should be used for the calculation.
result := DoEchelonMatTransformationBlockwise(
    A,
    rec(
        numberBlocksHeight := numberBlocks,
        numberBlocksWidth := numberBlocks
    )
);

# Of course you can specify an arbitrary combination of the arguments:
result := DoEchelonMatTransformationBlockwise(
    A,
    rec(
        galoisField := GF(q),
        numberBlocksHeight := numberBlocks,
        numberBlocksWidth := numberBlocks
    )
);
```

In general the result of `DoEchelonMatTransformationBlockwise` is a record
containing the `vectors`, `pivotrows`, `pivotcols`, `rank` and `heads`. If
`withTrafo := true` then it also contains `coeffs`, `relations` and
`transformation`.

## Benchmark Data

To get some benchmark data on
- how the algorithm performs in relation to the `Gauss` pkg
- wall and CPU time
- [lock contention](
https://en.wikipedia.org/wiki/Lock_%28computer_science%29#Granularity)

start HPC-GAP and execute
```
Read("read.g");
n := 4000;; numberBlocks := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
MeasureContention(numberBlocks, q, A);
```

The algorithm will print out information and statistics.

## Compatibility with old HPC-GAP version
Older versions of hpcgap still have working guards. Unfortunately, in those
versions you can't load this package via `LoadPackage`.
Let's assume that you can run your old HPC-GAP version via $HPCGAP-OLD
Instead `cd` into the root folder of the `GaussPar` package and run:
`$HPCGAP-OLD compatibility-for-old-hpcgap/read.g`

## Contact

Please submit bug reports, suggestions for improvements and patches via
the [issue tracker](https://github.com/lbfm-rwth/GaussPar/issues)
or via email to
[Sergio Siccha](mailto:siccha@mathematik.uni-kl.de)
or
[Jendrik Brachter](mailto:brachter@cs.uni-kl.de).

## License

`GaussPar` is free software you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.
See file gpl-2.0.txt for further information.
