#
# GaussParallel: Parallel gaussian algorithm for finite fields
#
# This file is a script which compiles the package manual.
#
if fail = LoadPackage("AutoDoc", "2018.02.14") then
    Error("AutoDoc version 2018.02.14 or newer is required.");
fi;
LoadPackage("GaussParallel");

scan_dirs := [
    "doc",
    "gap",
    "gap/benchmarking"
    ];

AutoDoc(rec(
    autodoc := rec( scan_dirs := scan_dirs ),
    gapdoc := rec(scan_dirs := scan_dirs),
    scaffold := true,
));
QUIT;
