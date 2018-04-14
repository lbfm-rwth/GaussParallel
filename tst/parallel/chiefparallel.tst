gap> Read("read_hpc.g");
gap> Read("main_full_par_trafo.g");
gap> n := 200;; numberChops := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> result := ChiefParallel(GF(q), A, numberChops, numberChops);;
gap> result_std := EchelonMatTransformation(A);;
gap> -1 * result.vectors = result_std.vectors;
true
gap> -1 * result.coeffs = result_std.coeffs;
true
