gap> Read("read.g");
gap> Read("main_seq_trafo.g");
gap> Read("echelon_form.g");
gap> dimension := 10;; rank := 5;; q := 5;; numberChops := 2;;
gap> rs := RandomSource(IsMersenneTwister);;
gap> echelon := echelonMat(dimension, dimension, rank, rs, GF(q));;
gap> shapeless := shapelessMat(echelon, dimension, dimension, rs, GF(q));;
gap> result := Chief(GF(q), shapeless, numberChops, numberChops);;
gap> result_std := EchelonMatTransformation(shapeless);;
gap> -1 * result.vectors = result_std.vectors;
true
gap> -1 * result.coeffs = result_std.coeffs;
true
gap> -Concatenation(result.coeffs, result.relations) * shapeless = echelon;
true
