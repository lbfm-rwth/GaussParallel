Chief := function( galoisField,mat,a,b )
    local   tmp,
            bs,
            A,
            B,
            C,
            D,
            E,
            K,
            M,
            R,
            X,
            i,
            j,
            k,
            k_,
            l,
            h;
    ##Preparation: Init and chopping the matrix mat
    C := ChopMatrix( galoisField,mat,a,b );
    A := [];
    B := [];
    D := [];
    E := [];
    K := [];
    M := [];
    R := [];
    X := [];
    for i in [ 1 .. a ] do
        A[i] := [];
        E[i] := [];
        K[i] := [];
        for k in [ 1 .. b ] do
            A[i][k] := rec( A:=[],M:=[],K:=[],rho:=[],E:=[],lambda:=[] );
            E[i][k] := rec( rho:=[],delta:=[],nr:=0 );
        od;
        for h in [ 1 .. a ] do
            K[i][h] := [];
        od;
    od;
    for i in [ 1 .. b ] do
        M[i] := [];
        for k in [ 1 .. a ] do
            M[i][k] := [];
        od;
    od;
    for k in [ 1 .. b ] do
        D[k] := rec( remnant:=[],bitstring:=[] );
        B[k] := [];
        R[k] := [];
        X[k] := [];
        for j in [ 1 .. b ] do
            B[k][j] := [];
            R[k][j] := [];
            X[k][j] := [];;
        od;
    od;
    ###############################
    ###############################

    ## Step 1 ##
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            tmp := ClearDown( galoisField,C[i][j],D[j],i );
            D[j] := tmp.D;
            A[i][j] := tmp.A;
            if j=1 then
                bs := rec( rho:=[],delta:=[],nr:=0 );
            else
                bs := E[i][j-1];
            fi;
            E[i][j] := Extend( A[i][j],bs,j );
            for k in [ j+1 .. b ] do
                tmp := UpdateRow(  galoisField,A[i][j],C[i][k],
                                    B[j][k],i );
                C[i][k] := tmp.C;
                B[j][k] := tmp.B;
            od;
            for h in [ 1 .. i ] do
                tmp := UpdateRowTrafo(  galoisField,A[i][j],K[i][h],
                                    M[j][h],E[h][j],i,h,j );
                K[i][h] := tmp.K;
                M[j][h] := tmp.M;
            od;
        od;    
    od;
    
    ## Step 2 ##
    for j in [ 1 .. b ] do
        for h in [ 1 .. a ] do
            M[j][h] := RowLengthen( galoisField,M[j][h],E[h][j],E[h][b] );            
        od;
    od;
    
    ## Step3 ##
    for k in [ 1 .. b ] do
        R[k][k] := D[k].remnant;
    od;
    for k_ in [ 1 .. b ] do
        k := b-k_+1;
        for j in [ 1 .. (k - 1) ] do
            tmp := CEX( galoisField,D[k].bitstring,B[j][k] );
            X[j][k] := tmp[1];
            if IsEmpty(X[j][k]) then continue; fi;
            R[j][k] := tmp[2];
            for l in [ k .. b ] do
                if not IsEmpty(R[k][l]) then
                    R[j][l] := R[j][l] + X[j][k]*R[k][l];
                fi;
            od;
            for h in [ 1 .. a ] do
                if not IsEmpty(M[k][h]) then
                    M[j][h] := M[j][h] + X[j][k]*M[k][h];
                fi;
            od;
        od;
    od;

    return rec( remnant := R,transformation := M );

end;

Chief_debug := function( galoisField,mat,a,b )
    local   tmp,
            bs,
            A,
            B,
            C,
            D,
            E,
            K,
            M,
            R,
            X,
            i,
            j,
            k,
            h;
    ##Preparation: Init and chopping the matrix mat
    C := ChopMatrix( galoisField,mat,a,b );
    A := [];
    B := [];
    D := [];
    E := [];
    K := [];
    M := [];
    R := [];
    X := [];
    for i in [ 1 .. a ] do
        A[i] := [];
        E[i] := [];
        K[i] := [];
        for k in [ 1 .. b ] do
            A[i][k] := rec( A:=[],M:=[],K:=[],rho:=[],E:=[],lambda:=[] );
            E[i][k] := rec( rho:=[],delta:=[],nr:=0 );
        od;
        for h in [ 1 .. a ] do
            K[i][h] := [];
        od;
    od;
    for i in [ 1 .. b ] do
        M[i] := [];
        for k in [ 1 .. a ] do
            M[i][k] := [];
        od;
    od;
    for k in [ 1 .. b ] do
        D[k] := rec( remnant:=[],bitstring:=[] );
        B[k] := [];
        R[k] := [];
        X[k] := [];
        for j in [ 1 .. b ] do
            B[k][j] := [];
            R[k][j] := [];
            X[k][j] := [];;
        od;
    od;
    ###############################
    ###############################

    ## Step 1 ##
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            Error( "BEFORE CLEARDOWN" );

            tmp := ClearDown( galoisField,C[i][j],D[j],i );
            D[j] := tmp.D;
            A[i][j] := tmp.A;
            if j=1 then
                bs := rec( rho:=[],delta:=[],nr:=0 );
            else
                bs := E[i][j-1];
            fi;
            E[i][j] := Extend( A[i][j],bs,j );
            Error( "AFTER CLEARDOWN" );
            for k in [ j+1 .. b ] do
                tmp := UpdateRow(  galoisField,A[i][j],C[i][k],
                                    B[j][k],i );
                C[i][k] := tmp.C;
                B[j][k] := tmp.B;
            od;
            Error( "AFTER UPDATEROW" );
            for h in [ 1 .. i ] do
                tmp := UpdateRowTrafo(  galoisField,A[i][j],K[i][h],
                                    M[j][h],E[h][j],i,h,j );
                K[i][h] := tmp.K;
                M[j][h] := tmp.M;
            od;
            Error( "AFTER UPDATETRAFO" );
        od;    
    od;

    return true;

end;
