gap> ReadPackage("GaussPar", "tst/testdata/matrices.g");;
gap> ReadPackage("GaussPar", "tst/testfunctions.g");;
gap> for i in [1..6] do result := GAUSS_TestSpecialMatrices(M[1], M_width[1], M_height[1], randomSource, GF(M_q[1]), M_numberChops[1], false); if not result then Print("Error: Special matrix number i"); fi; od;
gap> # No matrix
gap> echelon := M[7];;
gap> shapeless := echelon;;
gap> result := DoEchelonMatTransformationBlockwise(shapeless, GF(M_q[7]), true, M_numberChops[7], M_numberChops[7]);
Error, no method found! For debugging hints type ?Recovery from NoMethodFound
Error, no 1st choice method found for `DimensionsMat' on 1 arguments
