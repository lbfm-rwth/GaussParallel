# Install global functions of this package

InstallGlobalFunction( DoEchelonMatTransformationBlockwise,
    function ( mat, options )
    # options is a record that can optionally specify
    # galoisField, IsHPC(removed for now), numberBlocksHeight, numberBlocksWidth,
    # withTrafo and verify
        local galoisField, numberBlocksHeight, numberBlocksWidth,
            dim, recnames, withTrafo, verify;

        recnames := Set( RecNames( options ) );

        dim := DimensionsMat( mat );
        if ("numberBlocksHeight" in recnames) and ("numberBlocksWidth" in recnames) then
            numberBlocksHeight := options.numberBlocksHeight;
            numberBlocksWidth := options.numberBlocksWidth;
        else
            numberBlocksHeight := GAUSS_calculateBlocks( dim[1] );
            numberBlocksWidth := GAUSS_calculateBlocks( dim[2] );
        fi;

        if "withTrafo" in recnames then
            withTrafo := options.withTrafo;
        else
            withTrafo := true;
        fi;

        if "galoisField" in recnames then
            galoisField := options.galoisField;
        else
            galoisField := DefaultFieldOfMatrix( mat );
            if galoisField = fail then
                Error("Please specify the field of the matrix using the options parameter.");
            fi;
        fi;

        # verify only makes sense if withTrafo is true
        if "verify" in recnames then
            verify := options.verify;
        else
            verify := withTrafo;
        fi;
        if verify and not withTrafo then
            Error("can't verify the calculation without computing the ",
                  "transformation matrix");
        fi;

        Info(InfoGauss, 1, "The matrix is split into ", numberBlocksHeight,
            " blocks vertically and ", numberBlocksWidth, " horizontally.");

        return Chief( galoisField, mat, numberBlocksHeight, numberBlocksWidth, withTrafo, verify );
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
