##############################################################################
##
#W  teststandard.g                   Gauss Package                    Sergio Siccha
##
##

Read("read.g");
TestDirectory("tst/standard", rec(exitGAP := true));
FORCE_QUIT_GAP(1);
