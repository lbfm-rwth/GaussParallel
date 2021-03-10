# This file contains functions that are called by functions from "gap/tasks.g".
# In principle, these functions only work in a single thread-local region and
# don't need to know anything about other threads.
# These functions may only access read-only objects or objects from the
# executing thread's thread-local region. They may only emit or write into
# thread-local objects.

# REX: RowEXtract
# mat is a matrix or an empty list
# positionsBitstring is a list of 0s or 1s.
# Generates two new matrices or possibly empty lists up, down from the rows of
# mat.
# The rows of `up` and `down` are copies of the rows of mat.
# 0 means a row goes into `down`, 1 means `up`.
GAUSS_REX := function(galoisField, positionsBitstring, mat)
    local down, allCols, upIndices, downIndices, up;
    if IsEmpty(mat) then
        return [[], []];
    fi;
    upIndices := Positions(positionsBitstring, 1);
    downIndices := Positions(positionsBitstring, 0);
    if IsEmpty(upIndices) then
        down := MutableCopyMat(mat);
        ConvertToMatrixRepNC(down, galoisField);
        return [[], down];
    fi;
    if IsEmpty(downIndices) then
        up := MutableCopyMat(mat);
        ConvertToMatrixRepNC(up, galoisField);
        return [up, []];
    fi;
    allCols := [1..NrCols(mat)];
    up := ExtractSubMatrix(mat, upIndices, allCols);
    down := ExtractSubMatrix(mat, downIndices, allCols);
    ConvertToMatrixRepNC(up, galoisField);
    ConvertToMatrixRepNC(down, galoisField);
    return [up, down];
end;

# CEX: ColumnEXtract
# Does the same as REX but on the columns of mat.
GAUSS_CEX := function(galoisField, positionsBitstring, mat)
    local right, allRows, leftIndices, rightIndices, left,
        transposed;
    if IsEmpty(mat) then
        return [[], []];
    fi;
    leftIndices := Positions(positionsBitstring, 1);
    rightIndices := Positions(positionsBitstring, 0);
    if IsEmpty(leftIndices) then
        right := MutableCopyMat(mat);
        ConvertToMatrixRepNC(right, galoisField);
        return [[], right];
    fi;
    if IsEmpty(rightIndices) then
        left := MutableCopyMat(mat);
        ConvertToMatrixRepNC(left, galoisField);
        return [left, []];
    fi;
    allRows := [1..NrRows(mat)];
    left := ExtractSubMatrix(mat, allRows, leftIndices);
    right := ExtractSubMatrix(mat, allRows, rightIndices);
    ConvertToMatrixRepNC(left, galoisField);
    ConvertToMatrixRepNC(right, galoisField);
    return [left, right];
end;

# PVC: PiVotCombine
# s, t are bitstrings.
# s is a bitstring of pivot rows or columns
# t is a shorter bitstring of pivots that we want merge into s.
# The entries of t overwrite the 0-entries of s
# Returns
# - newBitstring: an updated copy of s.
# - u: a bitstring of length "number 1's in newBitstring" that indicates which
#      1's come from s or t respectively
GAUSS_PVC := function (s, t)
    local newBitstring, u, positionU, positionT, i;
    #Case s is empty
    if IsEmpty(s) then
        u := ListWithIdenticalEntries(Sum(t), 1);
        return [t, u];
    fi;
    if IsEmpty(t) then
        u := ListWithIdenticalEntries(Sum(s), 0);
        return [s, u];
    fi;
    #Case otherwise
    u := [];
    newBitstring := [];
    positionU := 1;
    positionT := 1;
    for i in [1..Length(s)] do
        if s[i] = 1 then
            newBitstring[i] := 1;
            u[positionU] := 0;
            positionU := positionU + 1;
        else
            if t[positionT] = 1 then
                newBitstring[i] := 1;
                u[positionU] := 1;
                positionU := positionU + 1;
            else
                newBitstring[i] := 0;
            fi;
            positionT := positionT + 1;
        fi;
    od;
    return [newBitstring, u];
end;

# RRF: RowRiFfle
# u is a list of 0s and 1s
# rows0, rows1 are matrices
# Number of rows of rows0 is number of 0s in u
# Number of rows of rows1 is number of 1s in u
# Generates a new matrix by combining the rows of row0 and row1 according to u
GAUSS_RRF := function(galoisField, rows0, rows1, u)
    local l, dim, index, index0, index1, new;
    if IsEmpty(rows0) then
        ConvertToMatrixRepNC(rows1, galoisField);
        return rows1;
    fi;
    if IsEmpty(rows1) then
        ConvertToMatrixRepNC(rows0, galoisField);
        return rows0;
    fi;
    l := Length(u);
    index0 := 1;
    index1 := 1;
    new := NullMat(l,NrCols(rows0));
    index := 0;
    while (index0 + index1 -1 <= l) do
        index := index0 + index1 -1;
        if u[index] = 0 then
            new[index] := ShallowCopy(rows0[index0]);
            index0 := index0 + 1;
        else
            new[index] := ShallowCopy(rows1[index1]);
            index1 := index1 + 1;
        fi;
    od;

    ConvertToMatrixRepNC(new, galoisField);
    return new;
end;

# CRZ: ColumnRiffleZero
# nr is a positive integer
# mat is a matrix, empty or with nr many rows
# u is a list of 0s and 1s
# Generates a new matrix by combining the columns of mat and zero-columns
#   according to u. For each 1 in u a zero-column is added to mat.
# Uses RRF on the Transposed input as a subroutine
GAUSS_CRZ := function(galoisField, mat, u, nr)
    local nullMat, numZero, tmp, sum;
    if IsEmpty(u) then return mat; fi;
    sum := Sum(u);
    if IsEmpty(mat) then
        return [];
    fi;
    nullMat := NullMat(sum, NrRows(mat), galoisField);
    ConvertToMatrixRepNC(nullMat, galoisField);
    return TransposedMat(GAUSS_RRF(galoisField, TransposedMat(mat), nullMat,
    u));
end;

# ADI: ADdIdentity
# mat is a matrix
# bitstring is a list of 0s and 1s
# Writes 1s into the columns of mat indexed by the positions of bitstring that
# contain a 1.
# Is used together with CRZ. Calling CRZ and ADI inserts columns of the
# identity matrix into mat.
GAUSS_ADI := function(galoisField, mat, bitstring)
    local copy, colIndices, one, i;
    if IsEmpty(mat) then
        copy := IdentityMatrix(Is8BitMatrixRep, galoisField, Length(bitstring));
        ConvertToMatrixRepNC(copy, galoisField);
        return copy;
    fi;
    copy := MutableCopyMat(mat);
    ConvertToMatrixRepNC(copy, galoisField);
    colIndices := Positions(bitstring, 1);
    one := One(galoisField);
    for i in [1 .. Length(colIndices)] do
        copy[i][colIndices[i]] := one;
    od;
    return copy;
end;

# MKR: MaKeRiffle
# bitstring, subBitstring are lists of 0s and 1s of the same length
# An entry subBitstring[i] can only be 1, if bitstring[i] is 1.
# Returns a new list newBitstring of length "number of 1s in bitstring"
#   indicating if subBitstring is 1 or 0 for the positions of bitstring equal
#   to 1
GAUSS_MKR := function(bitstring, subBitstring)
    local newBitstring, current, i, l;
    if IsEmpty(subBitstring) then
        return ListWithIdenticalEntries(Sum(bitstring), 1);
    fi;
    l := Length(bitstring);
    newBitstring := [];
    current := 1;
    for i in [1 .. l] do
        if bitstring[i] = 0 then
            continue;
        else
            if subBitstring[i] = 1 then
                newBitstring[current] := 0;
            else
                newBitstring[current] := 1;
            fi;
        current := current + 1;
        fi;
    od;
    return newBitstring;
end;


# ECH: Echelonize
# galoisField is a finite field
# H is a matrix
# Creates 3 new matrices: M, K, R and 2 new lists of 0s and 1s: rho, gamma
# R encodes the remnant of a row reduced echelon form (RREF) of H, i.e.
# ( -1  | R )
# (  0  | 0 )
# is the row reduced echelon form of H.
# The other return values encode the transformation to get there.
# For a definition of these objects refer to the paper. They satisfy the
# equality:
# ( M | 0 )  ( rho)                            ( -1  | R )
# ( K | 1 )  (!rho)  H  (gamma^T | !gamma^T) = (  0  | 0 )
#
# EchelonMatTransformation computes the RREF of H and a transformation to
# get there. We need to transform this output into M, K, R, rho, gamma.
# The return values of EchelonMatTransformation satisfy:
# (   coeffs  )       ( vectors )
# ( relations )  H  = (    0    )
#
# This means:
# ( M | 0 )  ( rho)   (   coeffs  )
# ( K | 1 )  (!rho) = ( relations )
#
# and
#
# ( vectors )                         ( -1  | R )
# (    0    )  (gamma^T | !gamma^T) = (  0  | 0 )
#
# (note that permutation matrices are orthogonal).
#
# We only need to reorder and throw away columns.
GAUSS_ECH := function(galoisField, H)
    local EMT, R, vectors, M, coeffs, K, relations, rho, gamma, nrPivots, Id,
    selectedRho, currentPivot, nonSelectedGamma, i;

    # If H is empty or zero matrix we return early!
    if IsEmpty(H) or IsZero(H) then
        return ListWithIdenticalEntries(5, []);
    fi;

    EMT := EchelonMatTransformation(H);
    # R remnant
    # We create this from vectors
    R := [];
    ConvertToMatrixRepNC(EMT.vectors, galoisField);
    vectors := TransposedMat(EMT.vectors);
    # M is the transformation to get the non-zero rows of the RREF
    # We create this from coeffs
    M := [];
    ConvertToMatrixRepNC(EMT.coeffs, galoisField);
    coeffs := TransposedMat(EMT.coeffs);
    # K is the transformation to get the zero rows of the RREF
    # We create this from relations
    K := [];
    ConvertToMatrixRepNC(EMT.relations, galoisField);
    relations := TransposedMat(EMT.relations);
    # rho is a row-select list. Is used to select the pivot, here non-zero, rows.
    rho := [];
    # gamma is a column-select list. Is used to select the pivot columns.
    gamma := [];
    # Our identity matrix has as many rows as R.
    nrPivots := NrCols(vectors);
    Id := IdentityMatrix(Is8BitMatrixRep, galoisField, nrPivots);
    ConvertToMatrixRepNC(Id, galoisField);
    # FIXME use heads?
    # (M|0) x P_rho = coeffs;
    for i in [1..NrRows(coeffs)] do
        if IsZero(coeffs[i]) then
            rho[i] := 0;
        else
            rho[i] := 1;
        fi;
    od;
    selectedRho := Positions(rho, 1);
    M := ExtractSubMatrix(coeffs, selectedRho, [1..NrCols(coeffs)]);
    if not IsEmpty(relations) then
        K := ExtractSubMatrix(relations, selectedRho, [1..NrCols(relations)]);
    fi;
    # Look for rows of the identity matrix
    currentPivot := 1;
    for i in [1..NrRows(vectors)] do
        if vectors[i] = Id[currentPivot] then
            gamma[i] := 1;
            currentPivot := currentPivot + 1;
        else
            gamma[i] := 0;
        fi;
        if currentPivot > nrPivots then
            break;
        fi;
    od;
    Append(gamma, ListWithIdenticalEntries(NrRows(vectors) - Length(gamma), 0));
    # R empty iff H full column rank
    nonSelectedGamma := Positions(gamma, 0);
    if IsEmpty(nonSelectedGamma) then
        R := [];
    else
        R := ExtractSubMatrix(vectors, nonSelectedGamma, [1..NrCols(vectors)]);
    fi;

    M := -TransposedMat(M);
    K := TransposedMat(K);
    R := -TransposedMat(R);
    return [M, K, R, rho, gamma];
 end;

# RowLengthen:
# mat is a matrix
# SubBitstring is a sub-bitstring of Bitstring as in the sense of GAUSS_MKR.
# The columns of mat correspond to 1s in SubBitstring
# The function creates a new matrix by
#   placing zero-columns in mat according to the positions which are
#   1 in Bitstring but not in SubBitstring (using GAUSS_CRZ)
# This is used to restore zero-columns in the transformation matrix that we
# didn't store explicitly during the algorithm.
GAUSS_RowLengthen := function(galoisField, mat, SubBitstring, Bitstring)
    local lambda;
    lambda := GAUSS_MKR(Bitstring.rho, SubBitstring.rho);
    return GAUSS_CRZ(galoisField, mat, lambda, SubBitstring.nr);
end;
