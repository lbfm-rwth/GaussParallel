# Install global functions of this package

InstallGlobalFunction( DoEchelonMatTransformationBlockwise,
    function ( mat, options )
    # options is a record that can optionally specify
    # galoisField, IsHPC, numberBlocksHeight, numberBlocksWidth
        local galoisField, IsHPC, numberBlocksHeight, numberBlocksWidth,
            dim, recnames;

        recnames := Set( RecNames( options ) );

        dim := DimensionsMat( mat );
        if ("numberBlocksHeight" in recnames) and ("numberBlocksWidth" in recnames) then
            numberBlocksHeight := options.numberBlocksHeight;
            numberBlocksWidth := options.numberBlocksWidth;
        else
            numberBlocksHeight := GAUSS_calculateBlocks( dim[1] );
            numberBlocksWidth := GAUSS_calculateBlocks( dim[2] );
        fi;
        Info(InfoGauss, 1, "The matrix is split into ", numberBlocksHeight,
            " blocks vertically and ", numberBlocksWidth, " horizontally.");
        if ((numberBlocksHeight = 1) or (numberBlocksWidth = 1)
            and (IsHPC = true)) then
            Info(InfoGauss, 1, "Warning: The size of the blocks is so small",
                " that the parallel version is unlikely to bring benefits",
                " in terms of runtime.");
        fi;

        if "IsHPC" in recnames then
            IsHPC := options.IsHPC;
        else
            if not IsHPCGAP then
                IsHPC := false;
            else
                IsHPC := true;
            fi;
        fi;

        if "galoisField" in recnames then
            galoisField := options.galoisField;
        else
            galoisField := DefaultFieldOfMatrix( mat );
            if galoisField = fail then
                return "Please specify the field of the matrix using the options parameter.";
            fi;
        fi;

        return Chief( galoisField, mat, numberBlocksHeight, numberBlocksWidth, IsHPC );
    end
);

InstallGlobalFunction( EchelonMatTransformationBlockwise,
    function ( mat )
        local result;
        result := DoEchelonMatTransformationBlockwise( mat, rec() );
        return rec(
            vectors := result.vectors,
            heads := result.heads,
            coeffs := result.coeffs,
            relations := result.relations
        );
    end
);

InstallGlobalFunction( EchelonMatBlockwise,
    function ( mat )
        local result;
        result := DoEchelonMatTransformationBlockwise( mat, rec() );
        return rec(
            vectors := result.vectors,
            heads := result.heads
        );
    end
);
