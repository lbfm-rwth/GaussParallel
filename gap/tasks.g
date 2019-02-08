# These functions are scheduled as tasks by the main routine. They need to make
# sure that they only write read-only objects into the "shared" atomic lists.

# E[i] contains information about the location of pivot-rows in block-row i
# Information in E[i,j] is based on information in E[i,j-1].
# Extend_destructive assembles the information in E[i,j] combining information
# from E[i][j-1] and the new generated infromation for block [i,j] which
# is contained in A[i,j]. 
# Writes into E[i,j]
GAUSS_Extend_destructive := function(A, E, i, j)
    local rhoE, nr, rhoA, res;
    if j = 1 then
        rhoE := MakeReadOnlyOrImmutableObj([]);
    else
        rhoE := E[i][j - 1].rho;
    fi;
    # Number of non-zero bits
    nr := Length(rhoE) - Sum(rhoE);
    rhoA := A[i][j].rho;
    res := GAUSS_PVC(rhoE, rhoA);
    res := rec(rho := res[1], delta := res[2], nr := nr);
    E[i][j] := MakeReadOnlyOrImmutableObj(res);
end;


# ClearDown_destructive is used to perform a local echelonization process in a
# singe block C[i,j] and to collect information about this process that is used
# by other tastks (e.g. location of pivot-rows/-cols, transformations, echelon
# form of block, ...)
# Writes into A[i,j] (information about pivot-rows and trafo) and 
#   D[j] (information about pivot-cols and the remnant in block-column j)
GAUSS_ClearDown_destructive := function(galoisField, C, D, A, i, j)
    local A_, M, K, bitstring, E, riffle, vectors_, vectors, H, tmp, ech;
    if IsEmpty(C[i][j]) then
        A[i][j] := MakeReadOnlyOrImmutableObj(
            rec(A := [], M := [], E := [], K := [], rho := [], lambda := [])
        );
        return;
    fi;
    if i = 1 then
        ech :=  GAUSS_ECH(galoisField, C[i][j]);
        tmp := GAUSS_PVC(MakeReadOnlyOrImmutableObj([]), ech[5]);

        A[i][j] := MakeReadOnlyOrImmutableObj(
            rec(A := [], M := ech[1], E := [], K := ech[2],
                rho := ech[4], lambda := tmp[2])
        );
        D[j] := MakeReadOnlyOrImmutableObj(
            rec(bitstring := ech[5], vectors := ech[3])
        );
        return;
    fi;
    tmp := GAUSS_CEX(galoisField, D[j].bitstring, C[i][j]);
    A_ := tmp[1];
    tmp := tmp[2];

    if IsEmpty(A_) or IsEmpty(D[j].vectors) then
        H := tmp;
    elif IsEmpty(tmp) then
        H := A_*D[j].vectors;
    else
        H := tmp + A_*D[j].vectors;
    fi;
    ech := GAUSS_ECH(galoisField, H);

    tmp := GAUSS_CEX(galoisField, ech[5], D[j].vectors);
    E := tmp[1];
    vectors_ := tmp[2];
    if not IsEmpty(ech[3]) and not IsEmpty(E) then
        vectors_ := vectors_ + E*ech[3];
    fi;
    tmp := GAUSS_PVC(D[j].bitstring, ech[5]);
    bitstring := tmp[1];
    riffle := tmp[2];
    vectors := GAUSS_RRF(galoisField, vectors_, ech[3], riffle);

    A[i][j] := MakeReadOnlyOrImmutableObj(
        rec(A := A_, M := ech[1], E := E, K := ech[2],
            rho := ech[4], lambda := riffle)
    );
    D[j] := MakeReadOnlyOrImmutableObj(
        rec(bitstring := bitstring, vectors := vectors)
    );
end;

# UpdateRow_destructive takes information from ClearDown(i,j) about 
# row transformations that are used to echelonize C[i,j] and propagates it to
# matrices to the right (i.e C[i,k] for a k > j).  The function basically
# copies these transformations and splits up C[i,k] into pivot-rows (going into
# B, not processed any further) and non-pivot-rows (staying in C) processed by
# other tasks later on.  
# Writes into C[i,j], B[j,k]
GAUSS_UpdateRow_destructive := function(galoisField, A, C, B, i, j, k)
    local tmp, B_, C_, S, V, W, X, Z;
    if IsEmpty(A[i][j].A) or IsEmpty(B[j][k]) then
        Z := C[i][k];
    elif IsEmpty(C[i][k]) then
        Z := A[i][j].A * B[j][k];
    else
        Z := C[i][k] + A[i][j].A * B[j][k];
    fi;

    tmp := GAUSS_REX(galoisField, A[i][j].rho, Z);
    V := tmp[1];
    W := tmp[2];
    X := [];
    if not IsEmpty(A[i][j].M) and not IsEmpty(V) then
        X := A[i][j].M * V;
    fi;
    if i > 1 then
        if IsEmpty(A[i][j].E) or IsEmpty(X) then
            S := B[j][k];
        elif IsEmpty(B[j][k]) then
            S := A[i][j].E * X;
        else
            S := B[j][k] + A[i][j].E * X;
        fi;
        B_ := GAUSS_RRF(galoisField, S, X, A[i][j].lambda);
    else
        B_ := X;
    fi;

    if IsEmpty(A[i][j].K) or IsEmpty(V)   then
        C_ := W;
    elif IsEmpty(W) then
        C_ := A[i][j].K * V;
    else
        C_ := W + A[i][j].K * V;
    fi;

    C[i][k] := MakeReadOnlyOrImmutableObj(C_);
    B[j][k] := MakeReadOnlyOrImmutableObj(B_);
end;

# Calls GAUSS_CEX but writes into R[j,k] and returns X
# PreClearUp and ClearUp_destructive peform a classical backwards-elimination
# process.
GAUSS_PreClearUp := function(R, galoisField, D, B, j, k)
    local tmp;
    tmp := GAUSS_CEX(galoisField, D[k].bitstring, B[j][k]);
    MakeReadOnlyOrImmutableObj(tmp[1]);
    MakeReadOnlyOrImmutableObj(tmp[2]);
    R[j][k] := tmp[2];
    return tmp[1];
end;

# Propagates the transformations used to echelonize C[i,j] to the transformation
# matrices M, K. Basically works like UpdateRow_destructive but uses implicit
# information about the shape of M and K.  See documentation of
# UpdateRow_destructive and refer to paper for changes.
# Writes into M[j,h] and K[i,h]
GAUSS_UpdateRowTrafo_destructive := function(galoisField, A, K, M, E, i, h, j)
    local tmp, K_, M_, S, V, W, X, Z;
    # for the or delta empty part, cf. paper: want to know if beta' in the
    # descr of GAUSS_UpdateRowTrafo is 0.
    if (IsEmpty(A[i][j].A) and IsEmpty(A[i][j].M)) or IsEmpty(E[h][j].delta) then
        return;
    fi;
    if j > 1 then
        K_ := GAUSS_CRZ(galoisField, K[i][h], E[h][j].delta, E[h][j].nr);
    fi;

    if (not h=i) and j > 1 then
        if IsEmpty(M[j][h]) or IsEmpty(A[i][j].A) then
            Z := K_;
        elif IsEmpty(K_) then
            Z := A[i][j].A * M[j][h];
        else
            Z := K_ + A[i][j].A * M[j][h];
        fi;
    elif (not h=i) then
        if IsEmpty(M[j][h]) then
            Z := Zero(galoisField) * A[i][j].A;
        else
            Z := A[i][j].A * M[j][h];
        fi;
    elif j>1 then
        Z := K_;
    fi;

    if not (j = 1 and h = i) then
        tmp := GAUSS_REX(galoisField, A[i][j].rho, Z);
        V := tmp[1];
        W := tmp[2];
    fi;

    if (not j = 1) and h = i  then
        V := GAUSS_ADI(galoisField, V, E[h][j].delta);
    fi;

    if not (j = 1 and h = i) then
        if IsEmpty(V) or IsEmpty(A[i][j].M) then
            X := A[i][j].M;
        else
            X := A[i][j].M * V;
        fi;
    else
        X := A[i][j].M;
    fi;

    if not h=i then
        if IsEmpty(X) or IsEmpty(A[i][j].E) then
            S := M[j][h];
        else
            S := M[j][h] + A[i][j].E * X;
        fi;
    elif not i=1 then
        if IsEmpty(X) or IsEmpty(A[i][j].E) then
            S := [];
        else
            S := A[i][j].E*X;
        fi;
    fi;

    if  not (h = i and i = 1) then
        M_ := GAUSS_RRF(galoisField, S, X, A[i][j].lambda);
    else
        M_ := X;
    fi;

    if  not (h = i and j = 1) then
        if IsEmpty(V) or IsEmpty(A[i][j].K) then
            K_ := W;
        elif IsEmpty(W) then
            K_ := A[i][j].K * V;
        else
            K_ := W + A[i][j].K * V;
        fi;
    else
        K_ := A[i][j].K;
    fi;

    K[i][h] := MakeReadOnlyOrImmutableObj(K_);
    M[j][h] := MakeReadOnlyOrImmutableObj(M_);
end;

# Writes into R[j,l]
GAUSS_ClearUp_destructive := function(R, X, j, k, l)
    if IsEmpty(R[k][l]) or IsEmpty(X) then return; fi;
    if IsEmpty(R[j][l]) then
        R[j][l] := MakeReadOnlyOrImmutableObj(X * R[k][l]);
    else
        R[j][l] := MakeReadOnlyOrImmutableObj(R[j][l] + X * R[k][l]);
    fi;
end;
