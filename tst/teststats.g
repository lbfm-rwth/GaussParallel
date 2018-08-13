##############################################################################
##
#W  teststats.g                   Gauss Package                    Emma Ahrens
##
##

# For now, test files should make sure to read the correct read*.g
Read("read.g");
TestDirectory("tst/stats", rec(exitGAP := true));
FORCE_QUIT_GAP(1);
