#
# GaussParallel: Parallel gaussian algorithm for finite fields
#
# Reading the declaration part of the package.
#
DeclareInfoClass("InfoGauss");

ReadPackage("GaussParallel", "gap/main.gd");
ReadPackage("GaussParallel", "gap/echelon-form.gd");
ReadPackage("GaussParallel", "gap/RREF.gd");
ReadPackage("GaussParallel", "gap/upstream.gd");

if IsHPCGAP then
    ReadPackage( "GaussParallel", "gap/benchmarking/measure_contention.gd");
fi;
