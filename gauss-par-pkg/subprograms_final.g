ChopMatrix := function( f,A,nrows,ncols )
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
    AA := [];

    for  i  in [ 1 .. nrows-1] do
        AA[i] := [];
        for j in [ 1 .. ncols-1 ] do
            AA[i][j] := ImmutableMatrix(f,A{[(i-1)*a+1 .. i*a]}{[(j-1)*b+1 .. j*b]});
        od;
    od;
    AA[nrows] := [];
    for i in [ 1 .. nrows-1 ] do
        AA[i][ncols] := ImmutableMatrix(f,A{[(i-1)*a+1 .. i*a]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]});
    od;
    for j in [ 1 .. ncols-1 ] do
        AA[nrows][j] := ImmutableMatrix(f,A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(j-1)*b+1 .. j*b]});
    od;
    AA[nrows][ncols] := ImmutableMatrix(f,A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]});    
    return AA;
end;

Extend := function( A,E,flag )
    local   tmp,
            rho,
            delta;

    if flag = 1 then
        tmp := PVC( [],A.rho  );
    else
        tmp := PVC( E.rho,A.rho );
    fi;
    return rec( rho := tmp[1],delta := tmp[2],nr := (Length(E.rho-Sum(E.rho))) );
end;

RowLengthen := function( galoisField,mat,Einter,Efin )
    local   lambda;
    lambda := MKR( Efin.rho,Einter.rho );
    return CRZ( galoisField,mat,lambda,Einter.nr );
end;

ClearDown := function( galoisField,C,D,i )
    local   A,
            M,
            K,
            bitstring,
            E,
            riffle,
            remnant_,
            remnant,
            H,
            tmp,
            ech;
    if IsEmpty(C) then
        return rec( A := rec(A:=[],M:=[],E:=[],K:=[],
                    rho:=[],lambda:=0*[1..Sum(D.bitstring)] ), D:=D );
    fi;
    if i = 1 then
        ech :=  ECH( galoisField,C );
        tmp := PVC( [],ech[5] );
        return rec( A := rec( A:=[],M:=ech[1],K:=ech[2],rho:=ech[4],E:=[],lambda:=tmp[2] )
        ,D:= rec(bitstring := ech[5],remnant := ech[3] ) );
    fi;
    tmp := CEX( galoisField,D.bitstring,C );
    A := tmp[1];
    tmp := tmp[2];
    
    if IsEmpty(A) or IsEmpty(D.remnant) then
        H := tmp;
    elif IsEmpty(tmp) then
        H := A*D.remnant;
    else
        H := tmp + A*D.remnant;
    fi;
    ech := ECH( galoisField,H );
    
    tmp := CEX( galoisField,ech[5],D.remnant );
    E := tmp[1];
    remnant_ := tmp[2];
    if not IsEmpty(ech[3]) and not IsEmpty(E) then
        tmp := tmp + E*ech[3];
    fi;
    tmp := PVC( D.bitstring,ech[5] );
    bitstring := tmp[1];
    riffle := tmp[2];
    remnant := RRF( galoisField,remnant_,ech[3],riffle );
 #   return rec( A := [ A,ech[1],ech[2],ech[4],E,ech[5] ], D:= [ bitstring,remnant ] );

    return rec( A := rec( A:=A,M:=ech[1],K:=ech[2],rho:=ech[4],E:=E,lambda:=riffle )
        ,D:= rec(bitstring := bitstring,remnant := remnant ) );
end;

UpdateRow := function( galoisField,A,C,B,i )
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

    tmp := REX( galoisField,A.rho,Z );
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
        B_ := RRF( galoisField,S,X,A.lambda );
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

UpdateRowTrafo := function( galoisField,A,K,M,E,i,h,j )
    local   tmp,
            K_,
            M_,
            S,
            V,
            W,
            X,
            Z;

    if j > 1 then
        K_ := CRZ( galoisField,K,E.delta,E.nr );
    fi;

    if ( not h=i ) and j > 1 then
        if IsEmpty(M) then
            Z := K_;
        else
            Z := K_ + A.A*M;
        fi;
    elif ( not h=i  ) then
        if IsEmpty(M) then
            Z := []; 
        else
            Z := A.A*M;
        fi;
    elif j>1 then
        Z := K_;
    fi;

    if not (j = 1 and h = i) then
        tmp := REX( galoisField,A.rho,Z );
        V := tmp[1];
        W := tmp[2];
    fi;

    if ( not j = 1 ) and h = i  then
        V := ADI( galoisField,V,E.delta );
    fi;

    if not (j = 1 and h = i) then
        X := A.M*V;
    else
        X := A.M;
    fi;

    if not h=i then
        S := M+A.E*X;
    elif not i=1 then
        S := A.E*X;
    fi;

    if  not ( h = i and i = 1 ) then
        M_ := RRF( galoisField,S,X,A.lambda );
    else
        M_ := X;
    fi;
    
    if  not ( h = i and j = 1 ) then
        K_ := W + A.K*V;
    else
        K_ := A.K;
    fi;

    return rec( K:=K_,M:=M_ );
end;
