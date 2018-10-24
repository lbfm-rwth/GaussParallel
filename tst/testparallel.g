##############################################################################
##
#W  testparallel.g                   Gauss Package                    Sergio Siccha
##
##

LoadPackage("GaussPar");
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/parallel"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
