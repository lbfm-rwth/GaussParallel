##############################################################################
##
#W  testextreme.g                   Gauss Package                    Sergio Siccha
##
##

LoadPackage("GaussPar");
SetInfoLevel(InfoGauss, 0);
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/extreme"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
