LoadPackage("GaussPar");
SetInfoLevel(InfoGauss, 0);

dirs := [
    DirectoriesPackageLibrary("GaussPar", "tst/parallel"),
    DirectoriesPackageLibrary("GaussPar", "tst/standard")
];
if IsHPCGAP then
    # timing suite
    Add(dirs, DirectoriesPackageLibrary("GaussPar", "tst/stats"));
fi;

TestDirectory(dirs, rec(exitGAP := true));

# if we ever get here, there was an error in how GAP handles tests
FORCE_QUIT_GAP(1);
