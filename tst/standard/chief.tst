gap> n := 200;; numberBlocks := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> result := DoEchelonMatTransformationBlockwise(A, rec( galoisField := GF(q), IsHPC := false, numberBlocksHeight := numberBlocks, numberBlocksWidth := numberBlocks ));;
gap> result_std := EchelonMatTransformation(A);;
gap> result.vectors = result_std.vectors;
true
gap> result.coeffs = result_std.coeffs;
true
