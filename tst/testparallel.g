##############################################################################
##
#W  testparallel.g                   Gauss Package                    Sergio Siccha
##
##

# For now, test files should make sure to read the correct read*.g
Read("read.g");
TestDirectory("tst/parallel", rec(exitGAP := true));
FORCE_QUIT_GAP(1);
