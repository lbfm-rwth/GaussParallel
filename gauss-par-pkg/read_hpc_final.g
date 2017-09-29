#LoadPackage("Gauss");
if not IsBound( EchelonMatTransformationDestructive ) then
    Read("./hpc/gauss-upwards.gd");
fi;
LoadPackage("IO");
Read("./hpc/gauss-upwards.gi");
Read("./utils_final.g");
Read("./subprograms_final.g");
Read("./main_seq_trafo_final.g");
#Read("./main_hpc_trafo.g");
#Read("./main_semi_hpc.g");
