# This file creates a matrix in strict row echelon form to be used for testing and a function that gets a matrix and returns another matrix with the same strict row echelon form.

echelonMat := function(height, width, rank, randomSeed, ring)
   local echelonMat, rightCorner, i, j;

    rightCorner := RandomMat(randomSeed, rank, width-rank, ring);
    echelonMat := NullMat(height, width);

    for i in [1 .. rank] do
        echelonMat[i][i] := 1;
    od;
    echelonMat := echelonMat * One(ring);

    for i in [1 .. rank] do
        for j in [1 .. (width-rank)] do
            echelonMat[i][rank+j] := rightCorner[i][j];
        od;
    od;

    return echelonMat;
end;

shapelessMat := function(mat, height, width, randomSeed, ring)
    return RandomInvertibleMat(randomSeed, height, ring) * mat;
end;
