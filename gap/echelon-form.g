# This file creates a matrix in strict row echelon form to be used for testing
# and a function that gets a matrix and returns another matrix with the same
# strict row echelon form.

InstallGlobalFunction( RandomEchelonMat,
function (height, width, rank, randomSeed, ring)
   local echelonMat, rightCorner, i, j;

    rightCorner := RandomMat(randomSeed, rank, width-rank, ring);
    echelonMat := NullMat(height, width, ring);

    for i in [1 .. rank] do
        echelonMat[i][i] := One(ring);
    od;
    echelonMat := echelonMat;

    for i in [1 .. rank] do
        for j in [1 .. (width-rank)] do
            echelonMat[i][rank+j] := rightCorner[i][j];
        od;
    od;

    return echelonMat;
end
);

GAUSS_RandomMatFromEchelonForm := function(echelon, height)
    local i, grp, left;

    # PseudoRandom takes too much time to initialize PseudoRandomSeed for big
    # values of height. So we call the initialization ourselves with parameters
    # that make it cheaper to compute.
    grp := GL(height, DefaultFieldOfMatrix(echelon));
    i := Length(GeneratorsOfGroup(grp));
    Group_InitPseudoRandom(grp, i+3, 20);
    left := PseudoRandom(grp);

    return left * echelon;
end;
