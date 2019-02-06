# We can trust the older versions of HPC-GAP to correctly place guards. Thus we
# use the old version to test whether we place objects in the right regions.
# We do not support the timing functions with older HPC-GAP versions.

IsHPCGAP := true;
MakeReadOnlyObj(MakeImmutable(IsHPCGAP));

# Content from init.g
DeclareInfoClass("InfoGauss");
SetInfoLevel(InfoGauss, 1);
Read("gap/main.gd");

# Content from read.g
# The gauss pkg doesn't work under old versions of HPCGAP. But we can take
# the functions we need from the following files we copied.
Read("compatibility-old-hpc/gauss-upwards.gd");
Read("compatibility-old-hpc/gauss-upwards.gi");

Read("gap/RREF.g");
Read("gap/dependencies.g");
Read("gap/thread-local.g");
Read("gap/utils.g");
Read("gap/tasks.g");

Read("gap/main.gi");
Read("gap/echelon-form.g");

# HACK
# We need to overwrite MakeReadOnlyOrImmutableObj
if IsHPCGAP then
     # In older hpcgap
     #      MakeReadOnly is recursive
     #      MakeReadOnlyObj is not recursive
     # In current hpcgap
     #      MakeReadOnlyObj is recursive
     #      MakeReadOnlySingleObj is not recursive
     MakeReadOnlyOrImmutableObj := MakeReadOnly;
else
     MakeReadOnlyOrImmutableObj := MakeImmutable;
fi;


Info(InfoGauss, 1, "<< The package \"GaussPar\" is still in alpha stage! >>");
Info(InfoGauss, 1, "<< See the README.md for some usage examples.      >>");
