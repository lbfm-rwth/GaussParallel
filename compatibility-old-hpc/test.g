# Taken from PR #149
GAUSS_BasicTest := function(IsHPC)
    local n, numberBlocks, q, A, result, result_std;

    n := 200;
    numberBlocks := 8;
    q := 5;
    A := RandomMat(n, n, GF(q));;

    result := DoEchelonMatTransformationBlockwise(A, rec( galoisField := GF(q), IsHPC := IsHPC, numberBlocksHeight := numberBlocks, numberBlocksWidth := numberBlocks ));
    result_std := EchelonMatTransformation(A);

    return (result.vectors = result_std.vectors)
        and (result.coeffs = result_std.coeffs);
end;

Read("compatibility-old-hpc/read.g");

success := true;
success := success and GAUSS_BasicTest(false);
success := success and GAUSS_BasicTest(true);

if not success then
    FORCE_QUIT_GAP(1);
fi;
