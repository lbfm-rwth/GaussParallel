LoadPackage("GaussPar");

if IsHPCGAP then
    # parallel
    TestDirectory(
        DirectoriesPackageLibrary("GaussPar", "tst/parallel"),
        rec(exitGAP := false)
    );
    # timing suite
    TestDirectory(
        DirectoriesPackageLibrary("GaussPar", "tst/stats"),
        rec(exitGAP := true)
    );
else
    # standard
    TestDirectory(
        DirectoriesPackageLibrary("GaussPar", "tst/standard"),
        rec(exitGAP := true)
    );
fi;

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
