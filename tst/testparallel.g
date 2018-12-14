##############################################################################
##
#W  testparallel.g                   Gauss Package                    Sergio Siccha
##
##

LoadPackage("GaussPar");
SetInfoLevel(InfoGauss, 0);
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/parallel"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
