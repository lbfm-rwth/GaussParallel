#LoadPackage("Gauss");
if not IsBound( EchelonMatTransformationDestructive ) then
    Read("./hpc/gauss-upwards.gd");
fi;
LoadPackage("IO");
Read("./hpc/gauss-upwards.gi");
Read("./utils.g");
Read("./subprograms.g");
Read("./main_trafo.g");
Read("./main_hpc_trafo.g");
Read("./main_semi_hpc.g");
