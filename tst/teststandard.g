##############################################################################
##
#W  testall.g                   Gauss Package                    Sergio Siccha
##
##

# For now, test files should make sure to read the correct read*.g
LoadPackage("gauss");
TestDirectory("tst/standard", rec(exitGAP := true));
FORCE_QUIT_GAP(1);
