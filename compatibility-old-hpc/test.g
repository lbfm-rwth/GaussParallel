Read("compatibility-old-hpc/read.g");
Read("tst/testfunctions.g");

failure := false;
for i in [1..5] do
    failure := not GAUSS_BasicTestEchelonMatTransformationBlockwise(
        180, 6, 6, 5, true
    );
    if failure then break; fi;
    failure := not GAUSS_BasicTestEchelonMatTransformationBlockwise(
        180, 6, 6, 5, false
    );
    if failure then break; fi;
    failure := not GAUSS_BasicTestEchelonMatTransformationBlockwise(
        240, 6, 6, 5, true
    );
    if failure then break; fi;
    failure := not GAUSS_BasicTestEchelonMatTransformationBlockwise(
        240, 6, 6, 5, false
    );
    if failure then break; fi;
od;

if failure then
    Print("Error!\n");
    FORCE_QUIT_GAP(1);
fi;

Print("Test successfull!\n");
QUIT_GAP(0);
