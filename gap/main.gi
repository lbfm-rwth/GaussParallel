# Install global functions of this package

InstallGlobalFunction( DoEchelonMatTransformationBlockwise,
    function ( mat, options... )
    # options is a list that can optionally specify
    # galoisField, IsHPC, numberBlocksHeight, numberBlocksWidth
        local galoisField, IsHPC, numberBlocksHeight, numberBlocksWidth,
            dim;

        dim := DimensionsMat( mat );
        if Size(options) < 4 then
            numberBlocksHeight := GAUSS_calculateBlocks( dim[1] );
            numberBlocksWidth := GAUSS_calculateBlocks( dim[2] );
        else
            numberBlocksHeight := options[3];
            numberBlocksWidth := options[4];
        fi;
        Info(InfoGauss, 1, "The matrix is split into ", numberBlocksHeight,
            " blocks vertically and ", numberBlocksWidth, " horizontally.");
        if ((numberBlocksHeight = 1) or (numberBlocksWidth = 1)
            and (IsHPC = true)) then
            Info(InfoGauss, 1, "Warning: The size of the blocks is so small",
                " that the parallel version is unlikely to bring benefits",
                " in terms of runtime.");
        fi;

        if Size(options) < 2 then
            if not IsHPCGAP then
                IsHPC := false;
            else
                IsHPC := true;
            fi;
        else
            IsHPC := options[2];
        fi;

        if Size(options) < 1 then
            galoisField := DefaultFieldOfMatrix( mat );
            if galoisField = fail then
                return "Please specify the field of the matrix using the options parameter.";
            fi;
        else
            galoisField := options[1];
        fi;

        return Chief( galoisField, mat, numberBlocksHeight, numberBlocksWidth, IsHPC );
    end
);

InstallGlobalFunction( EchelonMatTransformationBlockwise,
    function ( mat )
        local result;
        result := DoEchelonMatTransformationBlockwise( mat );
        return rec(
            vectors := result.vectors,
            heads := result.pivotrows, #TODO ersetze mit heads
            coeffs := result.coeffs,
            relations := result.relations
        );
    end
);

InstallGlobalFunction( EchelonMatBlockwise,
    function ( mat )
        local result;
        result := DoEchelonMatTransformationBlockwise( mat );
        return rec(
            vectors := result.vectors,
            heads := result.pivotrows, #TODO ersetze mit heads
        );
    end
);
