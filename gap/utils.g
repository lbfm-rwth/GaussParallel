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
    local   i, j, rrem, crem, rowCnt, rowOverhead, colOverhead, colCnt, AA, a,
    b;

    rrem := NrRows(A) mod nrows;
    crem := NrCols(A) mod ncols;
    a := ( NrRows(A) - rrem ) / nrows; 
    b := ( NrCols(A) - crem ) / ncols; 
    ## the alogirthm tries to chop the matrix A in equally sized submatrices
    ## the basic size of each block will be axb, the remainder of the above division is then
    ## spread among the frist rrem rows/cram cols respectively, giving each block 1 addiitonal row and/or column

    ## create a matrix AA of size 'nrows x ncols' which stores all submatrices
    AA := FixedAtomicList(nrows, 0);
    
    rowCnt := 0;
    rowOverhead := 1;
    for i in [ 1 .. nrows ] do
        colCnt := 0;
        colOverhead := 1;
        if  i > rrem then rowOverhead := 0; fi;
        AA[i] := FixedAtomicList(ncols, 0);
        for j in [ 1 .. ncols ] do
            if  j > crem then colOverhead := 0; fi;
            AA[i][j] := A{[rowCnt+1 .. rowCnt+a+rowOverhead]}{[colCnt+1 .. colCnt+b+colOverhead]};
            ConvertToMatrixRepNC(AA[i][j],f);
            MakeReadOnlyOrImmutableObj(AA[i][j]);
            colCnt := colCnt + b + colOverhead;
        od;
        rowCnt := rowCnt + a + rowOverhead;
    od;

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

GAUSS_GlueMorK := function( galoisField, blockMat, nrBlocksHeigth, nrBlocksWidth, mat, E, v, destinationRows, destinationCols, flag)
    local destination, a, b, X, rows, i, j, tmpR, tmpC, tmp, nullMat, w, M;
    M := blockMat;
    if   destinationRows = 0  or destinationCols = 0 then
        return [];
    fi;
    a := nrBlocksHeigth;
    b := nrBlocksWidth;
    
    destination := NullMat(destinationRows, destinationCols, galoisField);
    ConvertToMatrixRepNC(destination, galoisField);

    ## TODO The following info could be rank(blockCol[j]) ??
    rows := [];
    for j in [ 1 .. a ] do
        rows[j] := 0;
        for i in [ 1 .. b ] do
            if  not IsEmpty(M[j][i]) then
                rows[j] := NrRows(M[j][i]);
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

        # TODO NrRows(mat) = Length(v) ??;
        w := NrRows(mat)/b;
        for i in [ 1 .. b ] do
            if IsEmpty(M[j][i]) then
                if IsEmpty(E[i][Length(E[i])].rho) then
                    tmp := w;
                else
                    tmp := Length( E[i][Length(E[i])].rho );
                fi;
                tmpC := tmpC + tmp; continue;
                #M[j][i] := NullMat( rows[j],tmp,galoisField );
            else
                nullMat := NullMat(
                    Length(E[i][Length(E[i])].rho) - NrCols(M[j][i]),
                    NrRows(M[j][i]), galoisField
                );
                ConvertToMatrixRepNC(nullMat, galoisField); 
                M[j][i] := TransposedMat(
                    GAUSS_RRF(galoisField, nullMat, TransposedMat(M[j][i]),
                              E[i][Length(E[i])].rho)
                );
            fi;
            destination{[tmpR .. tmpR + NrRows(M[j][i])-1 ]}
                {[tmpC .. tmpC + NrCols(M[j][i])-1 ]}
                := M[j][i];
            tmpC := tmpC + NrCols(M[j][i]);
        od;
        tmpR := tmpR + rows[j];
    od;

    # FIXME: convert?
    ConvertToMatrixRepNC(destination, galoisField);
    if flag = 1 then
        X := IdentityMat( destinationRows,galoisField );
        ConvertToMatrixRepNC(X, galoisField);
        return TransposedMat( GAUSS_RRF( galoisField,X,
            TransposedMat( GAUSS_CEX( galoisField,v,destination )[1] ),v ) );
    fi;
    return destination;
end;

GAUSS_GlueFromBlocks := function( galoisField, blockMat, LocalSizeInfo, LocalColInfo, riffle, flag, hack )
    local i, j, a, b, nullMat, newRows, newCols, newMat, rowCnt, colCnt, tmp;

    a := Length(blockMat); # NrRows of blockMat
    b := Length(blockMat[1]); # NrCols of Blockmat
    
    # create a new matrix into which the blocks are glued 
    # compute its size first
    newRows := Sum(Concatenation(LocalSizeInfo));
    newCols := Length(Concatenation(LocalSizeInfo));
    if flag = 1 then
        newRows := newCols - newRows;
    fi;
    if  newRows = 0 or newCols = 0 then
        return [];
    fi;
    newMat := NullMat(newRows, newRows, galoisField);
    
    rowCnt := 1;
    for i in [ 1 .. a ] do
        colCnt := 1;
        if Sum(LocalColInfo[i])=0 then
            continue; 
        fi; 
        for j in [ 1 .. b ] do
            if  IsEmpty(blockMat[i][j]) then       
                if  IsEmpty(LocalSizeInfo[j]) then
                    tmp := 0;
                else
                    tmp := Sum(LocalSizeInfo[j]);
                fi;
                colCnt := colCnt + tmp;
                continue;
            else
                newMat{[rowCnt .. rowCnt + NrRows(blockMat[i][j])-1 ]}
                    {[colCnt .. colCnt + NrCols(blockMat[i][j])-1 ]}
                    := blockMat[i][j];
                colCnt := colCnt + NrCols(blockMat[i][j]);
            fi;
        od;
        rowCnt := rowCnt + Sum(LocalColInfo[i]);
    od;
    
    ConvertToMatrixRepNC(newMat, galoisField);
    if flag = 1 then
        tmp := IdentityMat(newRows, galoisField);
    else
        tmp := Length(riffle) - Sum(Concatenation(LocalSizeInfo));
        if  tmp = 0 then return newMat; fi;
        tmp := NullMat(tmp,newRows,galoisField);
    fi;
    ConvertToMatrixRep(tmp, galoisField);
    return TransposedMat( GAUSS_RRF(galoisField, tmp, TransposedMat(newMat), riffle)
               );
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
                 rows[i] := NrRows(R[i][j]);
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
                     tmpC := tmpC + NrCols(R[1][j]);
                 fi;
                 continue;
             fi;
             C{[tmpR .. tmpR + NrRows(R[i][j])-1 ]}
                {[tmpC .. tmpC + NrCols(R[i][j])-1 ]}
                := R[i][j];
            tmpC := tmpC + NrCols(R[i][j]);
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


GAUSS_WriteOutput := function( mat,a,b,ncols,nrows,galoisField,D,R,M,E,K,withTrafo )
    local v, w, rank, tmp, B, C, heads, i, j, GlueR, result, localBitstrings, localCols;

    ## Write output
    Info(InfoGauss, 2, "Write output");
    # Begin with row-select bitstring named v
    v := [];
    rank := 0;
    w := 0*[1..(NrRows(mat)/a) ];
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
        #B := GAUSS_GlueMorK( galoisField, M, b, a, mat, E, v, rank, Length(v), 0 );
        localBitstrings := ListWithIdenticalEntries(a,[]);
        localCols := ListWithIdenticalEntries(b,[]);
        for  j in [ 1 .. b ] do
            localCols[j] := D[j].bitstring;
        od;
        for i in [ 1 .. a ] do
            localBitstrings[i] := E[i][b].rho;
        od;
        B := GAUSS_GlueFromBlocks(galoisField, M, localBitstrings, localCols, v, 0, NrRows(mat)/b);  
    fi;

    GlueR := GAUSS_GlueR(rank, ncols, galoisField, nrows, D, R, a, b);
    C := GlueR.C;
    w := GlueR.w;

    if withTrafo then
        ## Glue the blocks of K
        D := GAUSS_GlueMorK(galoisField, K, a, a, mat, E, v, Length(v)-rank, Length(v), 1 );
    fi;
    
    heads := GAUSS_createHeads(v, w, NrCols(mat));
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
