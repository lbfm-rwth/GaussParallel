# Install global functions of this package

InstallGlobalFunction( DoEchelonMatTransformationBlockwise,
    function ( mat, options )
    # options is a record that can optionally specify
    # galoisField, IsHPC, numberChopsHeight, numberChopsWidth
        local galoisField, IsHPC, numberChopsHeight, numberChopsWidth,
            dim, recnames;

        recnames := Set( RecNames( options ) );

        dim := DimensionsMat( mat );
        if ("numberChopsHeight" in recnames) and ("numberChopsWidth" in recnames) then
            numberChopsHeight := options.numberChopsHeight;
            numberChopsWidth := options.numberChopsWidth;
        else
            numberChopsHeight := GAUSS_calculateChops( dim[1] );
            numberChopsWidth := GAUSS_calculateChops( dim[2] );
        fi;
        Info(InfoGauss, 1, "The matrix is split into ", numberChopsHeight,
            " chops vertically and ", numberChopsWidth, " horizontally.");
        if ((numberChopsHeight = 1) or (numberChopsWidth = 1)
            and (IsHPC = true)) then
            Info(InfoGauss, 1, "Warning: The size of the chops is so small",
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

        return Chief( galoisField, mat, numberChopsHeight, numberChopsWidth, IsHPC );
    end
);

InstallGlobalFunction( EchelonMatTransformationBlockwise,
    function ( mat )
        local result;
        result := DoEchelonMatTransformationBlockwise( mat );
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
        result := DoEchelonMatTransformationBlockwise( mat );
        return rec(
            vectors := result.vectors,
            heads := result.heads
        );
    end
);
