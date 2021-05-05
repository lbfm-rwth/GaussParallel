#
# GaussPar: Parallel gaussian algorithm for finite fields
#
# Reading the implementation part of the package.
#

# The gauss pkg doesn't work under old versions of HPCGAP. But we can take
# the functions we need from the following files we copied.
gauss := LoadPackage("GAUSS");
if gauss = fail then
    ReadPackage( "GaussPar", "compatibility-old-hpc/gauss-upwards.gd");
    ReadPackage( "GaussPar", "compatibility-old-hpc/gauss-upwards.gi");
fi;
LoadPackage("IO");

ReadPackage( "GaussPar", "gap/hpcgap-mockups.g");

ReadPackage( "GaussPar", "gap/RREF.g");
ReadPackage( "GaussPar", "gap/dependencies.g");
ReadPackage( "GaussPar", "gap/thread-local.g");
ReadPackage( "GaussPar", "gap/utils.g");
ReadPackage( "GaussPar", "gap/tasks.g");

ReadPackage( "GaussPar", "gap/main.gi");
ReadPackage( "GaussPar", "gap/echelon-form.g");

if IsHPCGAP then
    ReadPackage( "GaussPar", "gap/benchmarking/timing.g");
    ReadPackage( "GaussPar", "gap/benchmarking/measure_contention.gi");
fi;
