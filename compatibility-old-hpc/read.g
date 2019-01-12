# We can trust the older versions of HPC-GAP to correctly place guards. So we
# also use the old version to test our algorithm.
# We do not support the timing functions with older HPC-GAP versions.

IsHPCGAP := true;
MakeReadOnlyObj(MakeImmutable(IsHPCGAP));

DeclareInfoClass("InfoGauss");
SetInfoLevel(InfoGauss, 1);
Read("gap/main.gd");

# The gauss pkg doesn't work under old versions of HPCGAP. But we can take
# the functions we need from the following files we copied.
if true then
    Read("compatibility-old-hpc/gauss-upwards.gd");
    Read("compatibility-old-hpc/gauss-upwards.gi");
fi;

if not IsHPCGAP then
    Read("gap/overload-hpcgap-functions-in-gap.g");
fi;

Read("gap/subprograms.g");
Read("gap/dependencies.g");

# Compatibility hack for old HPC-GAP
ErrorNoReturn := Error;
Read("gap/main.g");
Read("gap/main.gi");
Read("gap/echelon_form.g");

Info(InfoGauss, 1, "<< The package \"GaussPar\" is still in alpha stage! >>");
Info(InfoGauss, 1, "<< See the README.md for some usage examples.      >>");
