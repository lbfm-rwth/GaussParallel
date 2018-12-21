GAUSS_TestSpecialMatrices := function(echelon, height, width, randomSource, galoisField, numberBlocks_height, numberBlocks_width, IsHPC)
    local shapeless, result, result_std;

    shapeless := GAUSS_RandomMatFromEchelonForm(echelon, width, height, randomSource, galoisField);
    result := DoEchelonMatTransformationBlockwise(
        shapeless,
        rec( galoisField := galoisField,
             IsHPC := IsHPC,
             numberBlocksHeight := numberBlocks_height,
             numberBlocksWidth := numberBlocks_width )
    );
    result_std := EchelonMatTransformation(shapeless);
    
    return (-1 * result.vectors = result_std.vectors)
        and (-1 * result.coeffs = result_std.coeffs)
        and (-Concatenation(result.coeffs, result.relations) * shapeless = echelon);
end;
