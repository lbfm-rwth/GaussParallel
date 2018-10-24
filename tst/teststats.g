##############################################################################
##
#W  teststats.g                   Gauss Package                    Emma Ahrens
##
##

LoadPackage("GaussPar");
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/stats"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
