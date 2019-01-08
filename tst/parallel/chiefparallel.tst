gap> ReadPackage("GaussPar", "tst/testfunctions.g");;
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(200, 8, 5, true);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(100, 10, 11, true);
true
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(300, 10, 17, true);
true
