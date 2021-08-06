##############################################################################
##
#W  teststandard.g                   Gauss Package                    Sergio Siccha
##
##

LoadPackage("GaussParallel");
SetInfoLevel(InfoGauss, 0);
TestDirectory(
    DirectoriesPackageLibrary("GaussParallel", "tst/standard"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
