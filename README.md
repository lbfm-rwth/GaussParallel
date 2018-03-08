# A Parallel Gaussian Algorithm
## How to use
Do something like this to run an example.
```
Read("read_hpc.g");
Read("main_full_par_trafo.g");
n := 4000;; chopSize := 8;; q := 5;;
A := RandomMat(n, n, GF(q));;
ChiefParallel(GF(q), A, chopSize, chopSize);;
```

To get some benchmark data on
- how the algorithm performs in relation to the `Gauss` pkg
- wall and CPU time
- lock contention

do `Read("measure_contention.g");`. The file will tell you how to proceed.
