gap> dimension := 200;; rank := 150;; q := 5;; numberBlocks := 8;;
gap> rs := RandomSource(IsMersenneTwister);;
gap> echelon := RandomEchelonMat(dimension, dimension, rank, rs, GF(q));;
gap> shapeless := GAUSS_shapelessMat(echelon, dimension, dimension, rs, GF(q));;
gap> result := DoEchelonMatTransformationBlockwise(shapeless, rec( galoisField := GF(q), IsHPC := true, numberBlocksHeight := numberBlocks, numberBlocksWidth := numberBlocks ));;
gap> result_std := EchelonMatTransformation(shapeless);;
gap> -1 * result.vectors = result_std.vectors;
true
gap> -1 * result.coeffs = result_std.coeffs;
true
gap> -Concatenation(result.coeffs, result.relations) * shapeless = echelon;
true
