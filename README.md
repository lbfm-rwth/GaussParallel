# How to use

Note that this is not an actual gap package but is supposed to be included in the GAP package Gauss. Therefore you can't load the code by executing LoadPackage(). Instead, there are two different ways to use GaussPar.

## The Sequential Implementation

You can use this by reading `read.g` and `main_seq_trafo.g`:
```
Read("read.g");
Read("main_seq_trafo.g");
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := Chief(GF(q), A, numberChops, numberChops);;
```

`result` is a record with the components `coeffs, vectors, relations, pivotrows, pivotcols, rank`. To try it out directly type
```
result.coeffs;
```
The result's components `coeffs, vectors, relations, pivotrows, pivotcols, rank` mimic the output of the (also sequential) `EchelonMatTransformation` of the Gauss package.
For further information see `?EchelonMatTransformation`.

## The Parallel Implementation

Use the parallel implementation in a similar way. Note that you need to start HPC-GAP!
```
Read("read_hpc.g");
Read("main_full_par_trafo.g");
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := ChiefParallel(GF(q), A, numberChops, numberChops);;
```

As in the sequential version the result is a record with the components `coeffs, vectors, relations, pivotrows, pivotcols, rank`.

## Benchmark Data

To get some benchmark data on
- how the algorithm performs in relation to the `Gauss` pkg
- wall and CPU time
- [lock contention](https://en.wikipedia.org/wiki/Lock_%28computer_science%29#Granularity)

do `Read("measure_contention.g");`. The file will tell you how to proceed.