# Used in test files

GAUSS_testMatrix := function(
        matrix, q, numberBlocksHeight, numberBlocksWidth, IsHPC)
    local result, result_std;

    result := DoEchelonMatTransformationBlockwise(matrix, rec( galoisField := GF(q), numberBlocksHeight := numberBlocksHeight, numberBlocksWidth := numberBlocksWidth ));;
    result_std := EchelonMatTransformation(matrix);;

    return (result.vectors = result_std.vectors)
        and (result.coeffs = result_std.coeffs)
        and (result.heads = result_std.heads);
end;

GAUSS_doubleTestMatrix := function(
        matrix, echelon, q, numberBlocksHeight, numberBlocksWidth, IsHPC)
    local result, result_std;

    result := DoEchelonMatTransformationBlockwise(matrix, rec( galoisField := GF(q), numberBlocksHeight := numberBlocksHeight, numberBlocksWidth := numberBlocksWidth ));;
    result_std := EchelonMatTransformation(matrix);;

    return (result.vectors = result_std.vectors)
        and (result.coeffs = result_std.coeffs)
        and (Concatenation(result.coeffs, result.relations) * matrix = echelon)
        and (result.heads = result_std.heads);
end;

GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank := function(
        dimension, rank, q, numberBlocksHeight, numberBlocksWidth, IsHPC)
    local echelon, shapeless, rs;
    rs := RandomSource(IsMersenneTwister);;

    echelon := RandomEchelonMat(dimension, dimension, rank, rs, GF(q));;
    shapeless := GAUSS_RandomMatFromEchelonForm(echelon, dimension);;

    return GAUSS_doubleTestMatrix(shapeless, echelon, q, numberBlocksHeight, numberBlocksWidth, IsHPC);
end;

GAUSS_TestEchelonMatTransformationBlockwiseWithGivenEchelonForm := function(
        echelon, height, width, randomSource, q, numberBlocksHeight,
        numberBlocksWidth, IsHPC)
    local shapeless, result, result_std;

    shapeless := GAUSS_RandomMatFromEchelonForm(echelon, width);
    return GAUSS_doubleTestMatrix(shapeless, echelon, q, numberBlocksHeight,
        numberBlocksWidth, IsHPC);
end;

GAUSS_BasicTestEchelonMatTransformationBlockwise := function(
        dimension, numberBlocksHeight, numberBlocksWidth, q, IsHPC)
    local matrix;

    matrix := RandomMat(dimension, dimension, GF(q));
    return GAUSS_testMatrix(matrix, q, numberBlocksHeight, numberBlocksWidth, IsHPC);
end;
