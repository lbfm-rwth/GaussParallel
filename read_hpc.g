# Loads all necessary functions into HPCGAP, assuming access to packages 
# "IO" and "GAUSS".

# Is this still true: My current versions of HPCGAP/ the GAUSS pkg are not 
# compatible, hence the repo at the moment
#            uses a rather obscure looking work-around..)


if not IsBound( EchelonMatTransformationDestructive ) then
    Read("./hpc/gauss-upwards.gd");
fi;
LoadPackage("IO");
Read("./hpc/gauss-upwards.gi");
Read("./utils.g");
Read("./subprograms.g");
Read("./dependencies_main.g");
Read("./main_seq_trafo.g");
Read("./main_full_par_trafo.g");
Read("./timing.g");
