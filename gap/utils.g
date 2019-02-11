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
    b, rowsList, colsList;

    rowsList := ListWithIdenticalEntries(nrows,0);
    colsList := ListWithIdenticalEntries(ncols,0);

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
        rowsList[i] := a+rowOverhead;
        AA[i] := FixedAtomicList(ncols, 0);
        for j in [ 1 .. ncols ] do
            if  j > crem then colOverhead := 0; fi;
            AA[i][j] := A{[rowCnt+1 .. rowCnt+rowsList[i]]}{[colCnt+1 .. colCnt+b+colOverhead]};
            ConvertToMatrixRepNC(AA[i][j],f);
            MakeReadOnlyOrImmutableObj(AA[i][j]);
            colsList[j] := b + colOverhead;

            colCnt := colCnt + b + colOverhead;
        od;
        rowCnt := rowCnt + a + rowOverhead;
    od;

    return rec(mat:=AA, rowsList:=rowsList, colsList:=colsList );
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


# This function takes one of the block-matrices M,K or R and writes the blocks
# into one big matrix. The cases for different block-matrices differ only in
# details, thus we have the flag to indicate we are work with R (flag = -1), 
# M (flag = 0) or K (flag = 1).
GAUSS_GlueFromBlocks := function( galoisField, blockMat, LocalSizeInfo, riffle,
flag)
    local i, j, a, b, nullMat, newRows, newCols, newMat, rowCnt, colCnt, tmp, rowInc;

    a := Length(blockMat); # NrRows of blockMat
    b := Length(blockMat[1]); # NrCols of Blockmat
    
    # create a new matrix into which the blocks are glued 
    # compute its size first
    if flag = 1 or flag = -1 then
        newCols := Sum(Concatenation(LocalSizeInfo));
        newRows := Length(riffle) - newCols;
    elif flag = 0 then
        newCols := Sum(Concatenation(LocalSizeInfo));
        newRows := newCols;
    fi;

    if newRows = 0 then
        return [];
    fi;
    if  newCols = 0 then
        if flag = 0 then
            return [];
        else
            return flag*IdentityMat(Length(riffle),galoisField);
        fi; 
    fi;
    newMat := NullMat(newRows, newCols, galoisField);
    
    rowCnt := 1;
    for i in [ 1 .. a ] do
        colCnt := 1;
        rowInc := 0;
        for j in [ 1 .. b ] do
            if  IsEmpty(blockMat[i][j]) then       
                if  IsEmpty(LocalSizeInfo[j]) then
                    tmp := 0;
                else
                    tmp := Sum(LocalSizeInfo[j]);
                fi;
                colCnt := colCnt + tmp;
            else
                rowInc := NrRows(blockMat[i][j]);
                newMat{[rowCnt .. rowCnt + NrRows(blockMat[i][j])-1 ]}
                    {[colCnt .. colCnt + NrCols(blockMat[i][j])-1 ]}
                    := blockMat[i][j];
                colCnt := colCnt + NrCols(blockMat[i][j]);
            fi;
        od;
        rowCnt := rowCnt + rowInc;
    od;
    
    ConvertToMatrixRepNC(newMat, galoisField);
    if not flag = 0 then
        tmp := flag*IdentityMat(newRows, galoisField);
    else
        tmp := Length(riffle) - Sum(Concatenation(LocalSizeInfo));
        if  tmp = 0 then return newMat; fi;
        tmp := NullMat(tmp,newRows,galoisField);
    fi; 
    ConvertToMatrixRep(tmp, galoisField);
    return TransposedMat( GAUSS_RRF(galoisField, tmp, TransposedMat(newMat), riffle)
               );
end;

GAUSS_WriteOutput := function(mat, a, b, nrColsPerBlockCol, nrRowsPerBlockRow,
                              galoisField, D, R, M, E, K, withTrafo)
    local v, rank, tmp, w, colBitstrings, rowBitstrings, B, C, heads, result,
    i, j;

    ## Write output
    Info(InfoGauss, 2, "Write output");
    # Begin with row-select bitstring named v
    v := [];
    rank := 0;

    for i in [ 1 .. a ] do
        if IsEmpty(E[i][b].rho) then
            tmp := ListWithIdenticalEntries(nrRowsPerBlockRow[i], 0);
        else
            tmp := E[i][b].rho;
        fi;
        v := Concatenation( v,tmp );
        rank := rank + Sum( tmp );
    od;

    w := [];
    for j in [ 1 .. b ] do
         if IsEmpty(D[j].bitstring) then
             D[j] := rec(remnant := [], bitstring := 0*[1..nrColsPerBlockCol[j]]);
             w := Concatenation(w,0*[1..nrColsPerBlockCol[j]]);
         else
             w := Concatenation(w,D[j].bitstring);
         fi;
    od;
    
    colBitstrings := ListWithIdenticalEntries(b,[]);
    for  j in [ 1 .. b ] do
        colBitstrings[j] := 1-D[j].bitstring;
    od;
    if  withTrafo then
        rowBitstrings := ListWithIdenticalEntries(a,[]);
        for i in [ 1 .. a ] do
            rowBitstrings[i] := E[i][b].rho;
        od;
        B := GAUSS_GlueFromBlocks(galoisField, M, rowBitstrings, v, 0);  
    fi;
    C := GAUSS_GlueFromBlocks(galoisField, R, colBitstrings, 1-w, -1);

    if withTrafo then
        D := GAUSS_GlueFromBlocks(galoisField, K, rowBitstrings, v, 1);
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
