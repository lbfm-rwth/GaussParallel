gap> dimension := 10;; rank := 5;; q := 5;; numberBlocks := 2;;
gap> rs := RandomSource(IsMersenneTwister);;
gap> echelon := RandomEchelonMat(dimension, dimension, rank, rs, GF(q));;
gap> shapeless := GAUSS_RandomMatFromEchelonForm(echelon, dimension);;
gap> result := DoEchelonMatTransformationBlockwise(shapeless, rec( galoisField := GF(q), IsHPC := false, numberBlocksHeight := numberBlocks, numberBlocksWidth := numberBlocks ));;
gap> result_std := EchelonMatTransformation(shapeless);;
gap> -1 * result.vectors = result_std.vectors;
true
gap> -1 * result.coeffs = result_std.coeffs;
true
gap> -Concatenation(result.coeffs, result.relations) * shapeless = echelon;
true
