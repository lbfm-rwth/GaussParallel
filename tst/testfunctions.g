GAUSS_TestSpecialMatrices := function(echelon, width, height, randomSource, galoisField, numberChops, IsHPC)
    local shapeless, result, result_std;

    shapeless := GAUSS_shapelessMat(echelon, width, height, randomSource, galoisField);
    result := DoEchelonMatTransformationBlockwise(shapeless, galoisField, IsHPC, numberChops, numberChops);
    result_std := EchelonMatTransformation(shapeless);
    
    return (-1 * result.vectors = result_std.vectors)
        and (-1 * result.coeffs = result_std.coeffs)
        and (-Concatenation(result.coeffs, result.relations) * shapeless = echelon);
end;
