#
# GaussPar: Parallel gaussian algorithm for finite fields
#
# Reading the declaration part of the package.
#
DeclareInfoClass("InfoGauss");
SetInfoLevel(InfoGauss, 1);

ReadPackage("GaussPar", "gap/main.gd");
ReadPackage("GaussPar", "gap/echelon-form.gd");
ReadPackage("GaussPar", "gap/RREF.gd");
ReadPackage("GaussPar", "gap/upstream.gd");

if IsHPCGAP then
    ReadPackage( "GaussPar", "gap/benchmarking/measure_contention.gd");
fi;
