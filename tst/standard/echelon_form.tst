gap> START_TEST("standard/echelon_form.tst");
gap> ReadPackage("GaussPar", "tst/testfunctions.g");;
gap> GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank(10, 2, 7, 5, 5, false);
true
gap> GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank(50, 30, 17, 2, 2, false);
true
gap> GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank(200, 157, 11, 8, 8, false);
true
gap> STOP_TEST("standard/echelon_form.tst");
