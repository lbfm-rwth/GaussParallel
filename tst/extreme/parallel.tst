gap> START_TEST("extreme/parallel.tst");
gap> ReadPackage("GaussParallel", "tst/testfunctions.g");;
gap> GAUSS_BasicTestEchelonMatTransformationBlockwise(4000, 35, 35, 7, true);
true
gap> STOP_TEST("extreme/parallel.tst");
