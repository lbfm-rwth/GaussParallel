gap> Read("read_hpc.g");
gap> Read("main_full_par_trafo.g");
gap> Read("measure_contention.g");
gap> n := 400;; numberChops := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> MeasureContention(numberChops, q, A, false);
