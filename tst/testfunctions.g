# Used in test files

GAUSS_testMatrix := function(matrix, q, numberBlocks_height, numberBlocks_width, IsHPC)
    local result, result_std;

    result := DoEchelonMatTransformationBlockwise(matrix, rec( galoisField := GF(q), IsHPC := IsHPC, numberBlocksHeight := numberBlocks_height, numberBlocksWidth := numberBlocks_width ));;
    result_std := EchelonMatTransformation(matrix);;

    return (result.vectors = result_std.vectors)
        and (result.coeffs = result_std.coeffs);
end;

GAUSS_doubleTestMatrix := function(matrix, echelon, q, numberBlocks_height, numberBlocks_width, IsHPC)
    local result, result_std;

    result := DoEchelonMatTransformationBlockwise(matrix, rec( galoisField := GF(q), IsHPC := IsHPC, numberBlocksHeight := numberBlocks_height, numberBlocksWidth := numberBlocks_width ));;
    result_std := EchelonMatTransformation(matrix);;

    return (result.vectors = result_std.vectors)
        and (result.coeffs = result_std.coeffs)
        and (Concatenation(result.coeffs, result.relations) * matrix = echelon);
end;

GAUSS_TestEchelonMatTransformationBlockwiseWithSetRank := function(dimension, rank, q, numberBlocks_height, numberBlocks_width, IsHPC)
    local echelon, shapeless, rs;
    rs := RandomSource(IsMersenneTwister);;

    echelon := RandomEchelonMat(dimension, dimension, rank, rs, GF(q));;
    shapeless := GAUSS_RandomMatFromEchelonForm(echelon, dimension);;

    return GAUSS_doubleTestMatrix(shapeless, echelon, q, numberBlocks_height, numberBlocks_width, IsHPC);
end;

GAUSS_TestEchelonMatTransformationBlockwiseWithGivenEchelonForm := function(echelon, height, width, randomSource, q, numberBlocks_height, numberBlocks_width, IsHPC)
    local shapeless, result, result_std;

    shapeless := GAUSS_RandomMatFromEchelonForm(echelon, width);
    return GAUSS_doubleTestMatrix(shapeless, echelon, q, numberBlocks_height, numberBlocks_width, IsHPC);
end;

GAUSS_BasicTestEchelonMatTransformationBlockwise := function(dimension, numberBlocks, q, IsHPC)
    local matrix;

    matrix := RandomMat(dimension, dimension, GF(q));
    return GAUSS_testMatrix(matrix, q, numberBlocks, numberBlocks, IsHPC);
end;
