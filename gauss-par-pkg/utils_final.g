REX := function( galoisField,positionsBitstring,mat )
    local   i,
            up,
            down,
            upCount,
            downCount,
            numrows,
            row;
    if IsEmpty ( mat ) then return [ [],[] ]; fi;
    numrows := Length( positionsBitstring );
    upCount := 1; 
    downCount := 1;
    up := []; 
    down := [];
    for i in [ 1 .. numrows ] do
        row := mat[ i ];
        if positionsBitstring[i] = 1 then
            up[ upCount ] := row; upCount := upCount + 1;
        else
            down[ downCount ] := row; downCount := downCount + 1;
        fi;
    od;
    ConvertToMatrixRepNC( up,galoisField );
    ConvertToMatrixRepNC( down,galoisField );

    return [ up,down ];
end;

CEX := function( galoisField,positionsBitstring,mat )
    local transposed;
    transposed := TransposedMat( mat );
    transposed := REX( galoisField,positionsBitstring,transposed );
    transposed[ 1 ] := TransposedMat( transposed[ 1 ] );
    transposed[ 2 ] := TransposedMat( transposed[ 2 ] );
    return transposed;
end;

PVC := function ( s,t )
    # We assume that positions of t correspond to zeroes in s
    local   newBitstring,
            u,
            positionU,
            positionT,
            i;
    u := []; 
    newBitstring := [];
    positionU := 1;
    positionT := 1;
    for i in [ 1..Length(s) ] do
        if s[ i ] = 1 then
            newBitstring[ i ] := 1;
            u[ positionU ] := 0;
            positionU := positionU + 1;
        else
            if t[ positionT ] = 1 then
                newBitstring[ i ] := 1;
                u[ positionU ] := 1;
                positionU := positionU + 1;
            else 
                newBitstring[ i ] := 0;
            fi;
            positionT := positionT + 1;
        fi;
    od;
    return [ newBitstring,u ]; 
end;
