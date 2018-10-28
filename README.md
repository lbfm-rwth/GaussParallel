# How to use

There are two different ways to use our package. You can either use the
sequential or the parallel version of the Gauss algorithm.
Both algorithms perform the Gauss algorithm on smaller chops of the
matrix and unite those chops later.

To use our package you should open a GAP session and type
```
LoadPackage("GaussPar");
```
with the directory GaussPar being somewhere in the paths searched through
by GAP.

Until now there are three different functions that you can use. The function `EchelonMatBlockwise` and `EchelonMatTransformationBlockwise` resemble `EchelonMat` and `EchelonMatTransformation` from the Gauss package. They need a matrix as input and return a record with different amounts of information.
The third function is called `DoEchelonMatTransformationBlockwise`. As input you can specify the matrix of course, the field of the matrix, the numbers of chops that are used and whether to use the sequential or the parallel version.

## EchelonMatBlockwise

The most basic usage of our algorithm would be the following:
```
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := EchelonMatBlockwise(A);;
```

`result` is a record with the components `vectors` and `heads`, whereas `vectors` is a list of the nonzero rows of the reduced echelon form of the matrix and `heads` is a list that of the pivot elements. If the i-th entry of `heads` is j, then we now that there is a pivot element at the position [j,i].

## EchelonMatTransformationBlockwise

To receive more information you should to this:
```
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
result := EchelonMatTransformationBlockwise(A);;
```

Now the record `result` does not contain only `vectors` and `heads`, but also `coeffs` and `relations`, whereas `coeffs` is the transformation matrix that the algorithm used to obtain the reduced row echelon form of the matrix and `relations` is the kernel of the matrix.

## DoEchelonMatTransformationBlockwise

There is another function that you can use when you want to specify a few things yourself:
```
n := 4000;; numberChops := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
# If you want to specify the field
result := DoEchelonMatTransformationBlockwise(A, GF(q));;

# If you want to specify, whether the sequential or parallel version
# is being used (note that you also have to specify the field in this case!)
# true means that you want to use the parallel version:
result := DoEchelonMatTransformationBlockwise(A, GF(q), true);

# And if you want to specify the number of chops, you have to specify
# the field and the version, too.
# The third argument stands for the number of chops that the height of
# the matrix is divided by and the fourth argument is the number of chops
# that the width of the matrix is divided by.
result := DoEchelonMatTransformationBlockwise(A, GF(q), true, numberChops, numberChops);
```

In those cases the record `result` contains every time `vectors`, `heads`, `coeffs` and `relations`. But also `pivotrows`, `pivotcols` and `rank` of the matrix.


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

The algorithm will out the statistics.

## Contact

Please submit bug reports, suggestions for improvements and patches via
the [issue tracker](https://github.com/lbfm-rwth/GaussPar/issues).

## License

`GaussPar` is free software you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.
