Read("read.g");

if IsHPCGAP then
    # parallel
    TestDirectory("tst/parallel", rec(exitGAP := false));
    # timing suite
    TestDirectory("tst/stats", rec(exitGAP := true));
else
    # standard
    TestDirectory("tst/standard", rec(exitGAP := true));
fi;

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
