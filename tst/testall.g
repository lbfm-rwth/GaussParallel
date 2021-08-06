LoadPackage("GaussParallel");
SetInfoLevel(InfoGauss, 0);

dirs := [
    DirectoriesPackageLibrary("GaussParallel", "tst/parallel"),
    DirectoriesPackageLibrary("GaussParallel", "tst/standard")
];
if IsHPCGAP then
    # timing suite
    Add(dirs, DirectoriesPackageLibrary("GaussParallel", "tst/benchmarking"));
fi;

TestDirectory(dirs, rec(exitGAP := true));

# if we ever get here, there was an error in how GAP handles tests
FORCE_QUIT_GAP(1);
