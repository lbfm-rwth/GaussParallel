# Install global functions of this package

InstallGlobalFunction( DoEchelonMatTransformationBlockwise,
    function ( mat, options... )
    # options is a list that can optionally specify
    # galoisField, IsHPC, numberChopsHeight, numberChopsWidth
        local galoisField, IsHPC, numberChopsHeight, numberChopsWidth,
            dim;

        dim := DimensionsMat( mat );
        if Size(options) < 4 then
            numberChopsHeight := GAUSS_calculateChops( dim[1] );
            numberChopsWidth := GAUSS_calculateChops( dim[2] );
        else
            numberChopsHeight := options[3];
            numberChopsWidth := options[4];
        fi;
        Info(InfoGauss, 1, "The matrix has a height of ", dim[1],
            " and a width of ", dim[2], " and is split into ", numberChopsHeight,
            " chops vertically and ", numberChopsWidth, " horizontally.");

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

        return Chief( galoisField, mat, numberChopsHeight, numberChopsWidth, IsHPC );
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
