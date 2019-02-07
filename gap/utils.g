# This function is temporary until it is not necessary anymore to find
# a number of blocks that divides the width or height.
# When different chop sizes are possible, we will presumably find the perfect
# chop size and divide in blocks of this size.
GAUSS_calculateBlocks := function( a )
    # a is either height or width of the matrix
    local i;
    i := 13;

    while (i > 1) do
        if (a mod i = 0) then
            return a/i;
        fi;
        i := i - 1;
    od;

    return 1;
end;

##############################################################################
# Larger high-level functions.
GAUSS_ChopMatrix := function( f,A,nrows,ncols )
    local   i,
            j,
            rrem,
            crem,
            AA,
            a,
            b;
    rrem := DimensionsMat(A)[1] mod nrows;
    crem := DimensionsMat(A)[2] mod ncols;
    a := ( DimensionsMat(A)[1] - rrem ) / nrows; 
    b := ( DimensionsMat(A)[2] - crem ) / ncols; 
    ## the alogirthm tries to chop the matrix A in equally sized submatrices
    ## create a matrix AA of size 'nrows x ncols' which stores all submatrices
    AA := FixedAtomicList(nrows, 0);
 
    ## these submatrices in AA have all equal dimensions 'a x b'
    for  i  in [ 1 .. nrows-1] do
        AA[i] := FixedAtomicList(ncols, 0);
        for j in [ 1 .. ncols-1 ] do
            AA[i][j] := A{[(i-1)*a+1 .. i*a]}{[(j-1)*b+1 .. j*b]};
            ConvertToMatrixRepNC(AA[i][j],f);
            MakeReadOnlyOrImmutableObj(AA[i][j]);
        od;
    od;

    ## to add the remaining submatrices we need to cut the submatrix dimensions if necessary
    AA[nrows] := FixedAtomicList(ncols, 0);
    for i in [ 1 .. nrows-1 ] do
        AA[i][ncols] := A{[(i-1)*a+1 .. i*a]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]};
        ConvertToMatrixRepNC(AA[i][ncols],f);
        MakeReadOnlyOrImmutableObj(AA[i][ncols]);
    od;
    for j in [ 1 .. ncols-1 ] do
        AA[nrows][j] := A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(j-1)*b+1 .. j*b]};
        ConvertToMatrixRepNC(AA[nrows][j],f);
        MakeReadOnlyOrImmutableObj(AA[nrows][j]);
    od;
    AA[nrows][ncols] := A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]};    
    ConvertToMatrixRepNC(AA[nrows][ncols],f);
    MakeReadOnlyOrImmutableObj(AA[nrows][ncols]);
    return AA;
end;

# Functions for Step 3
GAUSS_createHeads := function( pivotrows, pivotcols, width )
    # inputs: list that contains the row numbers of the pivot rows and
    # list that contains the column numbers of the pivot cols and
    # width of matrix
    local result, i, currentPivot;

    if not Length(pivotrows) = Length(pivotcols) then
        return [];
    fi;

    currentPivot := 0;
    result := ListWithIdenticalEntries( width, 0 );

    for i in [1 .. width] do
        if i in pivotcols then
            currentPivot := currentPivot + 1;
            result[i] := pivotrows[currentPivot];
        else
            result[i] := 0;
        fi;
    od;

    return result;
end;

GAUSS_GlueM := function(rank, v, galoisField, a, b, M, E, mat)
    local B, rows, i, j, tmpR, tmpC, tmp, nullMat, w;

    B := NullMat( rank,Length(v),galoisField );
    ConvertToMatrixRepNC(B, galoisField);
    rows := [];
    for j in [ 1 .. b ] do
        rows[j] := 0;
        for i in [ 1 .. a ] do
            if  not IsEmpty(M[j][i]) then
                rows[j] := DimensionsMat(M[j][i])[1];
                break;
            fi;
        od;
    od;

    tmpR := 1;
    for j in [ 1 .. b ] do
        if rows[j]=0 then
            continue;
        fi;
        tmpC := 1;
        w := DimensionsMat(mat)[1]/a;
        for i in [ 1 .. a ] do
            if IsEmpty(M[j][i]) then
                if IsEmpty(E[i][b].rho) then
                    tmp := w;
                else
                    tmp := Length( E[i][b].rho );
                fi;
                tmpC := tmpC + tmp; continue;
                #M[j][i] := NullMat( rows[j],tmp,galoisField );
            else
                # FIXME convert?
                nullMat := NullMat(Length(E[i][b].rho)
                                   - DimensionsMat(M[j][i])[2],
                                   DimensionsMat(M[j][i])[1],
                                   galoisField);
                M[j][i] := TransposedMat(
                    GAUSS_RRF(galoisField, nullMat, TransposedMat(M[j][i]),
                              E[i][b].rho)
                );
            fi;
            B{[tmpR .. tmpR + DimensionsMat(M[j][i])[1]-1 ]}
                {[tmpC .. tmpC + DimensionsMat(M[j][i])[2]-1 ]}
                := M[j][i];
            tmpC := tmpC + DimensionsMat(M[j][i])[2];
        od;
        tmpR := tmpR + rows[j];
    od;

    # FIXME: convert?
    return B;
end;

GAUSS_GlueR := function(rank, ncols, galoisField, nrows, D, R, a, b)
    local C, rows, w, i, j, tmpR, tmpC, idMat;

    # We can't convert C since `C{list1}{list2} := ..` is not valid for
    # compressed matrices. We need to use CopySubMatrix instead.
    C := NullMat( rank,ncols-rank,galoisField );
    rows := [];
    w := [];
    for i in [ 1 .. b ] do
         rows[i] := 0;
         if IsEmpty(D[i].bitstring) then
             w := Concatenation( w,0*[1..nrows[i]] );
         else
             w := Concatenation( w,D[i].bitstring );
         fi;
         for j in [ 1 .. b ] do
             if  not IsEmpty(R[i][j]) then
                 rows[i] := DimensionsMat(R[i][j])[1];
                 break;
             fi;
         od;
     od;

     tmpR := 1;
     for i in [ 1 .. b ] do
         if rows[i]=0 then
             continue;
         fi;
         tmpC := 1;
         for j in [ 1 .. b ] do
             if IsEmpty(R[i][j]) then
                 if not IsEmpty(D[j].bitstring) then
                     tmpC := tmpC + Sum( 1 - D[j].bitstring );
                 elif  not IsEmpty(R[1][j]) then
                     tmpC := tmpC + DimensionsMat(R[1][j])[2];
                 fi;
                 continue;
             fi;
             C{[tmpR .. tmpR + DimensionsMat(R[i][j])[1]-1 ]}
                {[tmpC .. tmpC + DimensionsMat(R[i][j])[2]-1 ]}
                := R[i][j];
            tmpC := tmpC + DimensionsMat(R[i][j])[2];
        od;
        tmpR := tmpR + rows[i];
    od;

    idMat := IdentityMat( rank,galoisField );
    ConvertToMatrixRepNC(idMat, galoisField);
    C := TransposedMat( GAUSS_RRF( galoisField, TransposedMat(C), -idMat, w ) );
    # FIXME Is it safe to call ConvertToMatrixRepNC(C, galoisField) here?
    # ConvertToMatrixRepNC(C, galoisField);
    return rec( C := C, w := w );
end;

GAUSS_GlueK := function(v, rank, galoisField, a, b, E, K, mat)
    local D, X, rows, i, j, tmpR, tmpC, tmp, nullMat;

    # We can't convert D and X since `D{list1}{list2} := ..` is not valid for
    # compressed matrices. We need to use CopySubMatrix instead.
    D := NullMat( Length(v)-rank,Length(v),galoisField );
    X := IdentityMat( Length(v)-rank,galoisField );
    rows := [];
    for j in [ 1 .. a ] do
        rows[j] := 0;
        for i in [ 1 .. a ] do
            if IsEmpty(K[j][i]) then
                rows[j] := Length(E[j][b].rho) - Sum(E[j][b].rho);
            else
                rows[j] := DimensionsMat(K[j][i])[1];
                break;
            fi;
        od;
    od;

    tmpR := 1;
    for j in [ 1 .. a ] do
        if rows[j]=0 then
            continue;
        fi;
        tmpC := 1;
        for i in [ 1 .. a ] do
            if IsEmpty(K[j][i]) then
                if IsEmpty(E[i][b].rho) then
                    #tmp := b;
                    tmp := DimensionsMat(mat)[1]/a;
                else
                    tmp := Length( E[i][b].rho );
                fi;
                tmpC := tmpC + tmp; continue;
                #K[j][i] := NullMat( rows[j],tmp,galoisField );
            else
                # FIXME convert?
                nullMat := NullMat(Length(E[i][b].rho)
                                   - DimensionsMat(K[j][i])[2],
                                   DimensionsMat(K[j][i])[1],
                                   galoisField);
                K[j][i] := TransposedMat(
                    GAUSS_RRF(galoisField, nullMat, TransposedMat(K[j][i]),
                              E[i][b].rho)
                );
            fi;
            D{[tmpR .. tmpR + DimensionsMat(K[j][i])[1]-1 ]}
                {[tmpC .. tmpC + DimensionsMat(K[j][i])[2]-1 ]}
                := K[j][i];
            tmpC := tmpC + DimensionsMat(K[j][i])[2];
        od;
        tmpR := tmpR + rows[j];
    od;

    # FIXME: Is it safe to call ConvertToMatrixRepNC(.., galoisField) on D and
    # X here?
    # D and X may contain empty lists as rows.
    # ConvertToMatrixRepNC(D, galoisField);
    # ConvertToMatrixRepNC(X, galoisField);
    return TransposedMat( GAUSS_RRF( galoisField,X,
        TransposedMat( GAUSS_CEX( galoisField,v,D )[1] ),v ) );
end;

GAUSS_WriteOutput := function( mat,a,b,ncols,nrows,galoisField,D,R,M,E,K,withTrafo )
    local v, w, rank, tmp, B, C, heads, i, j, GlueR, result;

    ## Write output
    Info(InfoGauss, 2, "Write output");
    # Begin with row-select bitstring named v
    v := [];
    rank := 0;
    w := 0*[1..(DimensionsMat(mat)[1]/a) ];
    for i in [ 1 .. a ] do
        if IsEmpty(E[i][b].rho) then
            tmp := w;
        else
            tmp := E[i][b].rho;
        fi;
        v := Concatenation( v,tmp );
        rank := rank + Sum( tmp );
    od;

    if  withTrafo then
        B := GAUSS_GlueM(rank, v, galoisField, a, b, M, E, mat);
    fi;

    GlueR := GAUSS_GlueR(rank, ncols, galoisField, nrows, D, R, a, b);
    C := GlueR.C;
    w := GlueR.w;

    if withTrafo then
        ## Glue the blocks of K
        D := GAUSS_GlueK(v, rank, galoisField, a, b, E, K, mat);
    fi;

    heads := GAUSS_createHeads(v, w, DimensionsMat(mat)[2]);
    result := rec(vectors := -C, pivotrows := v, pivotcols := w, rank := rank,
                  heads := heads);
    if withTrafo then
        result.coeffs := -B;
        result.relations := D;
        result.transformation := Concatenation(result.coeffs,
                                               result.relations);
    fi;
    
    return result;
end;
