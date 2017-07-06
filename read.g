#LoadPackage("Gauss");
if not IsBound( EchelonMatTransformationDestructive ) then
    Read("./hpcgap-version-gauss-upwards/gauss-upwards.gd");
fi;
LoadPackage("IO");
Read("./hpcgap-version-gauss-upwards/gauss-upwards.gi");
Read("./utils.g");
Read("./timing.g");
Read("./subprograms.g");
#Read("./step1_timed.g");
Read("./main_new.g");
Read("./main_new_parallel.g");
#Read("./step2.g");
#Read("./main.g");
