##############################################################################
##
#W  teststats.g                   Gauss Package                    Emma Ahrens
##
##

LoadPackage("GaussPar");
SetInfoLevel(InfoGauss, 0);
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/stats"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
