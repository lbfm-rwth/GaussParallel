##############################################################################
##
#W  teststats.g                   Gauss Package                    Emma Ahrens
##
##

LoadPackage("GaussParallel");
SetInfoLevel(InfoGauss, 0);
TestDirectory(
    DirectoriesPackageLibrary("GaussParallel", "tst/benchmarking"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
