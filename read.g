DeclareInfoClass("InfoGauss");
SetInfoLevel(InfoGauss, 0);

gauss := LoadPackage("GAUSS");
if gauss = fail then
    Read("./hpc/gauss-upwards.gd");
    Read("./hpc/gauss-upwards.gi");
fi;
LoadPackage("IO");

if not IsHPCGAP then
    Read("./tasks.g");
fi;

Read("./utils.g");
Read("./subprograms.g");
Read("./dependencies_main.g");

Read("./main.g");
Read("./timing.g");
Read("./echelon_form.g");

if IsHPCGAP then
    Read("./measure_contention.g");
    Read("./stats/timing.g");
fi;
