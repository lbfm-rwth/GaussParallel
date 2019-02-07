# This file contains functions that are called by functions from "gap/tasks.g".
# In principle, these functions only work in a single thread-local region and
# don't need to know anything about other threads.
# These functions may only access read-only objects or objects from the
# executing thread's thread-local region. They may only emit or write into
# thread-local objects.
GAUSS_REX := function(galoisField, positionsBitstring, mat)
    local i, up, down, upCount, downCount, numrows, row;
    if IsEmpty (mat) then
        return [[], []];
    fi;
    if IsEmpty (positionsBitstring) then
        return [[], mat];
    fi;
    numrows := Length(positionsBitstring);
    upCount := 1;
    downCount := 1;
    up := [];
    down := [];
    for i in [1 .. numrows] do
        row := ShallowCopy(mat[i]);
        if positionsBitstring[i] = 1 then
            up[upCount] := row; upCount := upCount + 1;
        else
            down[downCount] := row; downCount := downCount + 1;
        fi;
    od;
    ConvertToMatrixRepNC(up, galoisField);
    ConvertToMatrixRepNC(down, galoisField);
    return [up, down];
end;

GAUSS_CEX := function(galoisField, positionsBitstring, mat)
    local transposed;
    transposed := TransposedMat(mat);
    transposed := ShallowCopy(
        GAUSS_REX(galoisField, positionsBitstring, transposed)
    );
    transposed[1] := TransposedMat(transposed[1]);
    transposed[2] := TransposedMat(transposed[2]);
    return transposed;
end;

GAUSS_PVC := function (s, t)
    #We assume that positions of t correspond to zeroes in s
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
    new := [];
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

GAUSS_CRZ := function(galoisField, mat, u, nr)
    local nullMat, numZero, tmp, sum;
    if IsEmpty(u) then return mat; fi;
    sum := Sum(u);
    if IsEmpty(mat) then
        return [];
    fi;
    nullMat := NullMat(sum, DimensionsMat(mat)[1], galoisField);
    ConvertToMatrixRepNC(nullMat, galoisField);
    return TransposedMat(GAUSS_RRF(galoisField, TransposedMat(mat), nullMat,
    u));
end;

GAUSS_ADI := function(galoisField, mat, bitstring)
    local one, posOfNewOnes, i, copy, l;
    if IsEmpty(mat) then
        copy := IdentityMat(Length(bitstring), galoisField);
        ConvertToMatrixRepNC(copy, galoisField);
        return copy;
    fi;
    copy := MutableCopyMat(mat);
    l := Length(bitstring);
    one := One(galoisField);
    posOfNewOnes := [];
    for i in [1 .. l] do
        if bitstring[i] = 0 then
            continue;
        fi;
        Add(posOfNewOnes, i);
    od;
    l := Length(posOfNewOnes);
    for i in [1 .. l] do
        copy[i][posOfNewOnes[i]] := one;
    od;
    ConvertToMatrixRepNC(copy, galoisField);
    return copy;
end;

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

GAUSS_ECH := function(f, H)
    local sct, Mct, Kct, tct, Rct, EMT, m, k, M, K, R, S, N, r, s, t, i, ind,
    Id, one, zero, dims, dimId;

    if IsEmpty(H) then
        return ListWithIdenticalEntries(5, []);
    fi;

    EMT := EchelonMatTransformation(H);
    m := TransposedMat(EMT.coeffs);
    if IsEmpty(m) then
        return ListWithIdenticalEntries(5, []);
    fi;
    r := TransposedMat(EMT.vectors);
    k := TransposedMat(EMT.relations);
    s := [];
    t := [];
    R := [];
    M := [];
    K := [];
    one := One(f);
    zero := Zero(f);
    Id := IdentityMat(Length(EMT.vectors), f);
    ConvertToMatrixRepNC(Id, f);
    sct := 1;
    Mct := 1;
    Kct := 1;
    Rct := 1;
    tct := 1;

    ind := 1;
    dims := DimensionsMat(m);
    for i in [1..dims[1]] do
        if m[i] = zero*[1..dims[2]] then
            s[sct] := 0;
            sct := sct + 1;
        else
            M[Mct] := m[i];
            Mct := Mct + 1;
            if not IsEmpty(k) then
                K[Kct] := k[i];
                Kct := Kct + 1;
            fi;
            s[sct] := 1;
            sct := sct + 1;
        fi;
    od;

    ind := 1;
    dimId := DimensionsMat(Id);
    for i in [1..DimensionsMat(r)[1]] do
        if ind > dimId[1] then
            t[tct] := 0;
            tct := tct + 1;
            R[Rct] := r[i];
            Rct := Rct + 1;
            continue;
        fi;
        if r[i] = Id[ind] then
            t[tct] := 1;
            tct := tct + 1;
            ind := ind + 1;
        else
            t[tct] := 0;
            tct := tct + 1;
            R[Rct] := r[i];
            Rct := Rct + 1;
        fi;
    od;

    M := -TransposedMat(M);
    ConvertToMatrixRepNC(M, f);
    K := TransposedMat(K);
    ConvertToMatrixRepNC(K, f);
    R := -TransposedMat(R);
    ConvertToMatrixRepNC(R, f);
    return [M, K, R, s, t];
 end;

GAUSS_RowLengthen := function(galoisField, mat, Einter, Efin)
    local lambda;
    lambda := GAUSS_MKR(Efin.rho, Einter.rho);
    return GAUSS_CRZ(galoisField, mat, lambda, Einter.nr);
end;
