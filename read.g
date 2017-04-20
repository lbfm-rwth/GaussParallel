#LoadPackage("Gauss");
if not IsBound( EchelonMatTransformationDestructive ) then
    Read("./hpcgap-version-gauss-upwards/gauss-upwards.gd");
fi;
LoadPackage("IO");
Read("./hpcgap-version-gauss-upwards/gauss-upwards.gi");
Read("./utils.g");
Read("./timing.g");
Read("./main-seq.g");
Read("./step1_timed.g");
Read("./step2.g");
Read("./main.g");
