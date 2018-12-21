# This file contains the submodules of the parallel Gaussian algorithm.

##############################################################################
# Small low-level functions.
GAUSS_REX := function( galoisField,positionsBitstring,mat )
    local   i,
            up,
            down,
            upCount,
            downCount,
            numrows,
            row;
    if IsEmpty ( mat ) then
        return [ Immutable([]), Immutable([]) ];
    fi;
    if IsEmpty ( positionsBitstring ) then
        # FIXME: can this ever do something? It looks like mat always is immutable
        ConvertToMatrixRepNC( mat,galoisField );
        # FIXME: MakeReadOnlyObj doesn't move compressed matrices into the public
        # region. So mat always is immutable?
        MakeReadOnlyObj( mat );
        return [ Immutable([]),mat ];
    fi;
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
    # FIXME: Problem with MakeReadOnlyObj, see above
    MakeReadOnlyObj( up );
    ConvertToMatrixRepNC( down,galoisField );
    # FIXME: Problem with MakeReadOnlyObj, see above
    MakeReadOnlyObj( down );

    return [ up,down ];
end;

GAUSS_CEX := function( galoisField,positionsBitstring,mat )
    local transposed;
    transposed := TransposedMat( mat );
    transposed := ShallowCopy(
        GAUSS_REX(galoisField, positionsBitstring, transposed)
    );
    transposed[ 1 ] := TransposedMat( transposed[ 1 ] );
    transposed[ 2 ] := TransposedMat( transposed[ 2 ] );
    return MakeReadOnlyObj(MakeImmutable(transposed));
end;

GAUSS_PVC := function ( s,t )
    #We assume that positions of t correspond to zeroes in s
    local   newBitstring,
            u,
            positionU,
            positionT,
            i;
    #Case s is empty 
    if IsEmpty(s) then
        u := ListWithIdenticalEntries(Sum(t), 1);
        return MakeReadOnlyObj(MakeImmutable([ t,u ]));
    fi;
    if IsEmpty(t) then
        u := ListWithIdenticalEntries(Sum(s), 0);
        return MakeReadOnlyObj(MakeImmutable([ s,u ]));
    fi;
    #Case otherwise
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
    return MakeReadOnlyObj(MakeImmutable([ newBitstring,u ]));
end;

GAUSS_RRF := function( galoisField,rows0,rows1,u )
    local   l,
            dim,
            index,
            index0,
            index1,
            new;
    if IsEmpty(rows0) then 
        ConvertToMatrixRepNC( rows1,galoisField );
        return MakeReadOnlyObj(MakeImmutable(rows1));
    fi;
    if IsEmpty(rows1) then 
        ConvertToMatrixRepNC( rows0,galoisField );
        return MakeReadOnlyObj(MakeImmutable(rows0));
    fi;
    l := Length(u);
    index0 := 1;
    index1 := 1;
    new := [];
    index := 0;
    while ( index0 + index1 -1 <= l ) do
        index := index0 + index1 -1;
        if u[index] = 0 then
            new[index] := rows0[index0];
            index0 := index0 + 1;
        else 
            new[index] := rows1[index1];
            index1 := index1 + 1;
        fi;
    od;

    ConvertToMatrixRepNC( new,galoisField );
    return MakeReadOnlyObj(MakeImmutable(new));
end;

GAUSS_CRZ := function( galoisField,mat,u,nr ) 
#####################################
## !! Mind order of inputs to RRF !!
#####################################
    local   nullMat,
            numZero,
            tmp,
            sum;
    if IsEmpty(u) then return mat; fi;
    sum := Sum(u);
    if IsEmpty(mat) then
        if sum = 0 then
            return MakeReadOnlyObj(MakeImmutable([]));
        else
            return MakeReadOnlyObj(MakeImmutable(
                NullMat(nr, sum, galoisField)
            ));
        fi;
    fi;
    nullMat := MakeReadOnlyObj(MakeImmutable(
        NullMat(sum, DimensionsMat(mat)[1], galoisField)
    ));
    return MakeReadOnlyObj(MakeImmutable(
        TransposedMat(GAUSS_RRF(galoisField, TransposedMat(mat), nullMat, u))
    ));
end;

GAUSS_ADI := function( galoisField,mat,bitstring  )
    local   one,
            posOfNewOnes,
            i,
            copy,
            l;
    if IsEmpty(mat) then
        return IdentityMat( Length(bitstring),galoisField );
    else
        copy := MutableCopyMat( mat );
    fi;
    l := Length( bitstring );
    one := One( galoisField );
    posOfNewOnes := [];
    for i in [ 1 .. l ] do
        if bitstring[i] = 0 then
            continue;
        fi;
        Add( posOfNewOnes,i );
    od;
    l := Length( posOfNewOnes );
    for i in [ 1 .. l ] do
        copy[ i ][ posOfNewOnes[i] ] := one;
    od;
    ConvertToMatrixRepNC( copy,galoisField );
    return copy;
end;

GAUSS_MKR := function( bitstring,subBitstring )
    local   newBitstring,
            current,
            i,
            l;
    if IsEmpty(subBitstring) then
        return 0*[1..Sum(bitstring)]+1;
    fi;
    l := Length( bitstring  );
    newBitstring := [];
    current := 1;
    for i in [ 1 .. l ] do
        if bitstring[ i ] = 0 then
            continue;
        else
            if subBitstring[ i ] = 1 then
                newBitstring[ current ] := 0;
            else
                newBitstring[ current ] := 1;
            fi;
        current := current + 1;
        fi;
    od;
    return newBitstring;
end;

GAUSS_ECH := function( f,H )
    local   sct,
            Mct,
            Kct,
            tct,
            Rct,
            EMT,
            m,
            k,
            M,
            K,
            R,
            S,
            N,
            r,
            s,
            t,
            i,
            ind,
            Id,
            one,
            zero,
            dims,
            dimId;

    if IsEmpty(H) then
        return MakeReadOnlyObj(MakeImmutable(ListWithIdenticalEntries(5, [])));
    fi;
 
    EMT := EchelonMatTransformation( H );
    m := TransposedMat(EMT.coeffs);
    if IsEmpty(m) then
        return MakeReadOnlyObj(MakeImmutable(ListWithIdenticalEntries(5, [])));
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
    Id := IdentityMat( Length(EMT.vectors),f );
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
    ConvertToMatrixRepNC( M,f );
    K := TransposedMat(K);
    ConvertToMatrixRepNC( K,f );
    R := -TransposedMat(R);
    ConvertToMatrixRepNC( R,f );
    return MakeReadOnlyObj(MakeImmutable([M,K,R,s,t]));
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
            MakeReadOnlyObj(MakeImmutable(AA[i][j]));
        od;
    od;

    ## to add the remaining submatrices we need to cut the submatrix dimensions if necessary
    AA[nrows] := FixedAtomicList(ncols, 0);
    for i in [ 1 .. nrows-1 ] do
        AA[i][ncols] := A{[(i-1)*a+1 .. i*a]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]};
        ConvertToMatrixRepNC(AA[i][ncols],f);
        MakeReadOnlyObj(MakeImmutable(AA[i][ncols]));
    od;
    for j in [ 1 .. ncols-1 ] do
        AA[nrows][j] := A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(j-1)*b+1 .. j*b]};
        ConvertToMatrixRepNC(AA[nrows][j],f);
        MakeReadOnlyObj(MakeImmutable(AA[nrows][j]));
    od;
    AA[nrows][ncols] := A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]};    
    ConvertToMatrixRepNC(AA[nrows][ncols],f);
    MakeReadOnlyObj(MakeImmutable(AA[nrows][ncols]));
    return AA;
end;

GAUSS_Extend := function( A,E,flag )
    local   tmp,
            rho,
            delta;

    if flag = 1 then
        tmp := GAUSS_PVC( MakeReadOnlyObj(MakeImmutable([])),A.rho  );
    else
        tmp := GAUSS_PVC( E.rho,A.rho );
    fi;
    return rec( rho := tmp[1],delta := tmp[2],
    nr := (Length(E.rho-Sum(E.rho))) );
end;

GAUSS_RowLengthen := function( galoisField,mat,Einter,Efin )
    local   lambda;
    lambda := GAUSS_MKR( Efin.rho,Einter.rho );
    return GAUSS_CRZ( galoisField,mat,lambda,Einter.nr );
end;

GAUSS_ClearDown := function( galoisField,C,D,i )
    local   A,
            M,
            K,
            bitstring,
            E,
            riffle,
            vectors_,
            vectors,
            H,
            tmp,
            ech;

    if IsEmpty(C) then
        return rec( A := rec(A:=[],M:=[],E:=[],K:=[],
                    rho:=[],lambda:=[] ), D:=D );
    fi;
    if i = 1 then
        ech :=  GAUSS_ECH( galoisField,C );
        tmp := GAUSS_PVC( MakeReadOnlyObj(MakeImmutable([])),ech[5] );
        return rec( A := rec( A:=[],M:=ech[1],K:=ech[2],rho:=ech[4],E:=[],lambda:=tmp[2] )
        ,D:= rec(bitstring := ech[5],vectors := ech[3] ) );
    fi;
    tmp := GAUSS_CEX( galoisField,D.bitstring,C );
    A := tmp[1];
    tmp := tmp[2];
    
    if IsEmpty(A) or IsEmpty(D.vectors) then
        H := tmp;
    elif IsEmpty(tmp) then
        H := A*D.vectors;
    else
        H := tmp + A*D.vectors;
    fi;
    ech := GAUSS_ECH( galoisField,H );
    
    tmp := GAUSS_CEX( galoisField,ech[5],D.vectors );
    E := tmp[1];
    vectors_ := tmp[2];
    if not IsEmpty(ech[3]) and not IsEmpty(E) then
        vectors_ := vectors_ + E*ech[3];
    fi;
    tmp := GAUSS_PVC( D.bitstring,ech[5] );
    bitstring := tmp[1];
    riffle := tmp[2];
    vectors := GAUSS_RRF( galoisField,vectors_,ech[3],riffle );

    return rec( A := rec( A:=A,M:=ech[1],K:=ech[2],rho:=ech[4],E:=E,lambda:=riffle )
        ,D:= rec(bitstring := bitstring,vectors := vectors ) );
end;

GAUSS_ClearDown_destructive := function( galoisField,C,D,A,i,j )
    local   A_,
            M,
            K,
            bitstring,
            E,
            riffle,
            vectors_,
            vectors,
            H,
            tmp,
            ech;

    if IsEmpty(C[i][j]) then
        A[i][j] := MakeReadOnlyObj(MakeImmutable(
            rec(A := [], M := [], E := [], K := [], rho := [], lambda := [])
        ));
        return;
    fi;
    if i = 1 then
        ech :=  GAUSS_ECH( galoisField,C[i][j] );
        tmp := GAUSS_PVC(MakeReadOnlyObj(MakeImmutable([])), ech[5]);

        A[i][j] := MakeReadOnlyObj(MakeImmutable(
            rec(A := [], M := ech[1], E := [], K := ech[2],
                rho := ech[4], lambda := tmp[2])
        ));
        D[j] := MakeReadOnlyObj(MakeImmutable(
            rec(bitstring := ech[5], vectors := ech[3])
        ));
        return;
    fi;
    tmp := GAUSS_CEX( galoisField, D[j].bitstring, C[i][j] );
    A_ := tmp[1];
    tmp := tmp[2];
    
    if IsEmpty(A_) or IsEmpty(D[j].vectors) then
        H := tmp;
    elif IsEmpty(tmp) then
        H := A_*D[j].vectors;
    else
        H := tmp + A_*D[j].vectors;
    fi;
    MakeReadOnlyObj(MakeImmutable(H));
    ech := GAUSS_ECH( galoisField,H );
    
    tmp := GAUSS_CEX( galoisField,ech[5],D[j].vectors );
    E := tmp[1];
    vectors_ := tmp[2];
    if not IsEmpty(ech[3]) and not IsEmpty(E) then
        vectors_ := vectors_ + E*ech[3];
    fi;
    MakeReadOnlyObj(MakeImmutable(vectors_));
    tmp := GAUSS_PVC( D[j].bitstring,ech[5] );
    bitstring := tmp[1];
    riffle := tmp[2];
    vectors := GAUSS_RRF( galoisField,vectors_,ech[3],riffle );

    A[i][j] := MakeReadOnlyObj(MakeImmutable(
        rec(A := A_, M := ech[1], E := E, K := ech[2],
            rho := ech[4], lambda := riffle)
    ));
    D[j] := MakeReadOnlyObj(MakeImmutable(
        rec(bitstring := bitstring, vectors := vectors)
    ));
end;

GAUSS_UpdateRow := function( galoisField,A,C,B,i )
    local   tmp,
            B_,
            C_,
            S,
            V,
            W,
            X,
            Z;
    if IsEmpty(A.A) or IsEmpty(B) then
        Z := C;
    elif IsEmpty(C) then
        Z := A.A*B;
    else
        Z := C + A.A*B;
    fi;

    tmp := GAUSS_REX( galoisField,A.rho,Z );
    V := tmp[1];
    W := tmp[2];
    X := [];
    if not IsEmpty(A.M) and not IsEmpty(V) then
        X := A.M*V;
    fi;
    if i > 1 then
        if IsEmpty(A.E) or IsEmpty(X) then
            S := B;
        elif IsEmpty(B) then
            S := A.E*X;
        else 
            S := B + A.E*X;
        fi;
        B_ := GAUSS_RRF( galoisField,S,X,A.lambda );
    else
        B_ := X;
    fi; 

    if IsEmpty(A.K) or IsEmpty(V)   then
        C_ := W;
    elif IsEmpty(W) then
        C_ := A.K*V;
    else
        C_ := W + A.K*V;
    fi;

    return rec( C := C_,B := B_ );
end;

GAUSS_UpdateRow_destructive := function( galoisField,A,C,B,i,j,k )
    local   tmp,
            B_,
            C_,
            S,
            V,
            W,
            X,
            Z;
    if IsEmpty(A[i][j].A) or IsEmpty(B[j][k]) then
        Z := C[i][k];
    elif IsEmpty(C[i][k]) then
        Z := A[i][j].A*B[j][k];
    else
        Z := C[i][k] + A[i][j].A*B[j][k];
    fi;

    tmp := GAUSS_REX( galoisField,A[i][j].rho,Z );
    V := tmp[1];
    W := tmp[2];
    X := [];
    if not IsEmpty(A[i][j].M) and not IsEmpty(V) then
        X := A[i][j].M*V;
    fi;
    if i > 1 then
        if IsEmpty(A[i][j].E) or IsEmpty(X) then
            S := B[j][k];
        elif IsEmpty(B[j][k]) then
            S := A[i][j].E*X;
        else 
            S := B[j][k] + A[i][j].E*X;
        fi;
        B_ := GAUSS_RRF( galoisField,S,X,A[i][j].lambda );
    else
        B_ := X;
    fi; 

    if IsEmpty(A[i][j].K) or IsEmpty(V)   then
        C_ := W;
    elif IsEmpty(W) then
        C_ := A[i][j].K*V;
    else
        C_ := W + A[i][j].K*V;
    fi;

    C[i][k] := C_;
    B[j][k] := B_;
end;

GAUSS_UpdateRowTrafo := function( galoisField,A,K,M,E,i,h,j )
    local   tmp,
            K_,
            M_,
            S,
            V,
            W,
            X,
            Z;

    if (IsEmpty(A.A) and IsEmpty(A.M)) or IsEmpty(E.delta) then
        return rec( K:=K,M:=M );
    fi;  #### for the or delta empty part, cf. paper: want to know if beta' in the descr of GAUSS_UpdateRowTrafo is 0..
    if j > 1 then
        K_ := GAUSS_CRZ( galoisField,K,E.delta,E.nr );
    fi;

    if ( not h=i ) and j > 1 then
        if IsEmpty(M) or IsEmpty(A.A) then
            Z := K_;
        else
            Z := K_ + A.A*M;
        fi;
    elif ( not h=i  ) then
        if IsEmpty(M) then
            Z := Zero(galoisField)*A.A; 
        else
            Z := A.A*M;
        fi;
    elif j>1 then
        Z := K_;
    fi;

    if not (j = 1 and h = i) then
        tmp := GAUSS_REX( galoisField,A.rho,Z );
        V := tmp[1];
        W := tmp[2];
    fi;

    if ( not j = 1 ) and h = i  then
        V := GAUSS_ADI( galoisField,V,E.delta );
    fi;

    if not (j = 1 and h = i) then
        if IsEmpty(V) or IsEmpty(A.M) then
            X := A.M;
        else
            X := A.M*V;
        fi;
    else
        X := A.M;
    fi;

    if not h=i then
        if IsEmpty(X) or IsEmpty(A.E) then
            S := M;
        else
            S := M+A.E*X;
        fi;
    elif not i=1 then
        if IsEmpty(X) or IsEmpty(A.E) then
            S := [];
        else
            S := A.E*X;
        fi;
    fi;

    if  not ( h = i and i = 1 ) then
        M_ := GAUSS_RRF( galoisField,S,X,A.lambda );
    else
        M_ := X;
    fi;
    
    if  not ( h = i and j = 1 ) then
        if IsEmpty(V) or IsEmpty(A.K) then
            K_ := W;
        else
            K_ := W + A.K*V;
        fi;
    else
        K_ := A.K;
    fi;

    return rec( K:=K_,M:=M_ );
end;

# Calls GAUSS_CEX but writes into R and returns X
GAUSS_PreClearUp := function( R,galoisField,D,B,j,k )
    local tmp;
    tmp := GAUSS_CEX( galoisField, D[k].bitstring, B[j][k] );
    MakeReadOnlyObj(tmp[1]);
    MakeReadOnlyObj(tmp[2]);
    R[j][k] := tmp[2];
    return tmp[1];
end;

GAUSS_UpdateRowTrafo_destructive := function( galoisField,A,K,M,E,i,h,j )
    local   tmp,
            K_,
            M_,
            S,
            V,
            W,
            X,
            Z;

    if (IsEmpty(A[i][j].A) and IsEmpty(A[i][j].M)) or IsEmpty(E[h][j].delta) then
        return;
    fi;  #### for the or delta empty part, cf. paper: want to know if beta' in the descr of GAUSS_UpdateRowTrafo is 0..
    if j > 1 then
        K_ := GAUSS_CRZ( galoisField,K[i][h],E[h][j].delta,E[h][j].nr );
    fi;

    if ( not h=i ) and j > 1 then
        if IsEmpty(M[j][h]) or IsEmpty(A[i][j].A) then
            Z := K_;
        else
            Z := K_ + A[i][j].A*M[j][h];
        fi;
    elif ( not h=i  ) then
        if IsEmpty(M[j][h]) then
            Z := Zero(galoisField)*A[i][j].A; 
        else
            Z := A[i][j].A*M[j][h];
        fi;
    elif j>1 then
        Z := K_;
    fi;

    if not (j = 1 and h = i) then
        tmp := GAUSS_REX( galoisField,A[i][j].rho,Z );
        V := tmp[1];
        W := tmp[2];
    fi;

    if ( not j = 1 ) and h = i  then
        V := GAUSS_ADI( galoisField,V,E[h][j].delta );
    fi;

    if not (j = 1 and h = i) then
        if IsEmpty(V) or IsEmpty(A[i][j].M) then
            X := A[i][j].M;
        else
            X := A[i][j].M*V;
        fi;
    else
        X := A[i][j].M;
    fi;

    if not h=i then
        if IsEmpty(X) or IsEmpty(A[i][j].E) then
            S := M[j][h];
        else
            S := M[j][h]+A[i][j].E*X;
        fi;
    elif not i=1 then
        if IsEmpty(X) or IsEmpty(A[i][j].E) then
            S := [];
        else
            S := A[i][j].E*X;
        fi;
    fi;

    if  not ( h = i and i = 1 ) then
        M_ := GAUSS_RRF( galoisField,S,X,A[i][j].lambda );
    else
        M_ := X;
    fi;
    
    if  not ( h = i and j = 1 ) then
        if IsEmpty(V) or IsEmpty(A[i][j].K) then
            K_ := W;
        else
            K_ := W + A[i][j].K*V;
        fi;
    else
        K_ := A[i][j].K;
    fi;

    K[i][h]:=K_;
    M[j][h]:=M_;
end;

GAUSS_ClearUp := function( R,X,R_ )
    if IsEmpty(R_) or IsEmpty(X) then return R; fi;
    if IsEmpty(R) then
        return X*R_;
    else
        return R + X*R_;
    fi;
end;

# Writes into R
GAUSS_ClearUp_destructive := function( R,X,j,k,l )
    local threadLocalResult;
    # We use threadLocalResult to avoid a race condition:
    # maybe another task tries to access R[j][i] while it is not in the
    # read-only region.
    if IsEmpty(R[k][l]) or IsEmpty(X) then return; fi;
    if IsEmpty(R[j][l]) then
        threadLocalResult := X*R[k][l];
    else
        threadLocalResult := R[j][l] + X*R[k][l];
    fi;
    MakeReadOnlyObj(threadLocalResult);
    R[j][l] := threadLocalResult;
end;
