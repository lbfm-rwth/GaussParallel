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


## Output
At the moment the output looks different from what the `Gauss` pkg returns.
The correspondence of record component names is
(`Gauss` name = `GaussPar` name)
- vectors = remnant
- coeffs = transformation

---
# Old info
The following information need to be updated and should go - in the form of comments - into the files
they are about.

## File contents
### read_hpc.g
Loads all necessary functions into HPCGAP, assuming access to packages "IO" and "GAUSS".

Is this still true: My current versions of HPCGAP/ the GAUSS pkg are not compatible, hence the repo at the moment
             uses a rather obscure looking work-around..)

### main_full_par_trafo.g
Contains a version of the elimination alg. for HPCGAP computing RREF and a transformation, running completely in parallel.

### utils.g
Collection of small basic functions used in subfunctions of the algorithm

### subfunctions.g
Collection of larger subfunctions used in the Gaussian elimination alg.

## File contents of unused / unnecessary (?) files
### main_seq_trafo.g
Contains a version of the elimination alg. for standard GAP computing RREF and a transformation

### main_semi_par_trafo.g
Contains a version of the elimination alg. for HPCGAP computing RREF and a transformation, where the second step of the
                        algorithm runs in parallel ( using HPCGAP's task arch. )

### main_par_trafo.g
Contains a version of the elimination alg. for HPCGAP computing RREF and a transformation, where the first step of the
                        algorithm runs in parallel ( using HPCGAP's task arch. )

### read.g
Loads all necessary functions into standard GAP, assuming access to packages "IO" and "GAUSS".

