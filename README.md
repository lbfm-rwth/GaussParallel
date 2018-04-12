# How to use

First you need to note that this is not a actual gap package but is supposed to be included in the GAP package Gauss. That's why you can't load the code by executing LoadPackage(). Instead, there are two different ways to use to start using GaussPar.

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
For further information regarding the components `coeffs, vectors, relations, pivotrows, pivotcols, rank` and the original (sequential) Gauss package type e.g. `?EchelonMatTransformation`.

## The Parallel Implementation

Use the parallel implementation in a similar way:
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
- lock contention

do `Read("measure_contention.g");`. The file will tell you how to proceed.
