gap> n := 200;; numberChops := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> result := DoEchelonMatTransformationBlockwise(A, rec( galoisField := GF(q), IsHPC := true, numberChopsHeight := numberChops, numberChopsWidth := numberChops ));;
gap> result_std := EchelonMatTransformation(A);;
gap> -1 * result.vectors = result_std.vectors;
true
gap> -1 * result.coeffs = result_std.coeffs;
true
