LoadPackage("GaussPar");
SetInfoLevel(InfoGauss, 0);

if IsHPCGAP then
    # timing suite
    TestDirectory(
        DirectoriesPackageLibrary("GaussPar", "tst/stats"),
        rec(exitGAP := false)
    );
fi;

# parallel
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/parallel"),
    rec(exitGAP := false)
);

# standard
TestDirectory(
    DirectoriesPackageLibrary("GaussPar", "tst/standard"),
    rec(exitGAP := true)
);

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
