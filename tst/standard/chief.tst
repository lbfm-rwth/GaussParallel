gap> START_TEST("standard/chief.tst");
gap> ReadPackage("GaussPar", "tst/testfunctions.g");;
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(200, 8, 8, 5, false);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(100, 10, 10, 11, false);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(300, 10, 10, 17, false);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(180, 2, 3, 17, false);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(180, 3, 2, 17, false);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(180, 1, 3, 17, false);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(180, 2, 1, 17, false);
true
gap> STOP_TEST("standard/chief.tst");
