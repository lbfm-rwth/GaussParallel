gap> START_TEST("parallel/echelon_form.tst");
gap> ReadPackage("GaussPar", "tst/testfunctions.g");;
gap> GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank(10, 2, 7, 5, 5, true);
true
gap> GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank(50, 30, 17, 2, 2, true);
true
gap> GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank(200, 157, 11, 8, 8, true);
true
gap> STOP_TEST("parallel/echelon_form.tst");
