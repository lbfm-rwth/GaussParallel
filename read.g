#
# GaussParallel: Parallel gaussian algorithm for finite fields
#
# Reading the implementation part of the package.
#

# The gauss pkg doesn't work under old versions of HPCGAP. But we can take
# the functions we need from the following files we copied.
gauss := LoadPackage("GAUSS");
if gauss = fail then
    ReadPackage( "GaussParallel", "compatibility-old-hpc/gauss-upwards.gd");
    ReadPackage( "GaussParallel", "compatibility-old-hpc/gauss-upwards.gi");
fi;
LoadPackage("IO");

ReadPackage( "GaussParallel", "gap/hpcgap-mockups.g");

ReadPackage( "GaussParallel", "gap/RREF.g");
ReadPackage( "GaussParallel", "gap/dependencies.g");
ReadPackage( "GaussParallel", "gap/thread-local.g");
ReadPackage( "GaussParallel", "gap/utils.g");
ReadPackage( "GaussParallel", "gap/tasks.g");

ReadPackage( "GaussParallel", "gap/main.gi");
ReadPackage( "GaussParallel", "gap/echelon-form.g");

if IsHPCGAP then
    ReadPackage( "GaussParallel", "gap/benchmarking/timing.g");
    ReadPackage( "GaussParallel", "gap/benchmarking/measure_contention.gi");
fi;
