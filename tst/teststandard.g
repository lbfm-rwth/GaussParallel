##############################################################################
##
#W  teststandard.g                   Gauss Package                    Sergio Siccha
##
##

# For now, test files should make sure to read the correct read*.g
Read("read.g");
TestDirectory("tst/standard", rec(exitGAP := true));
FORCE_QUIT_GAP(1);
