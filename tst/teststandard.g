##############################################################################
##
#W  teststandard.g                   Gauss Package                    Sergio Siccha
##
##

LoadPackage("GaussPar");
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/standard"),
    rec(exitGAP := true)
);
FORCE_QUIT_GAP(1);
