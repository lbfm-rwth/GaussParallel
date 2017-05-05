GaussParallel := function( Inp,a,b ) #Chop inputmatrix Inp into (a)x(b) matrix
    local C,D,B,A,E,F,M,K,X,R, tmp,tmpR,tmpC,f,i,j,k,h,v,w,rank,rows,ncols;
    C := ChoppedMatrix( Inp,a,b );w := []; v :=[];R := [];A := [];B := [];D := [];E := [];F := [];M := [];K := [];X := [];
    f := DefaultFieldOfMatrix( Inp );	
    ncols := DimensionsMat( Inp )[2];
    Inp := MutableCopyMat(C); 
    # Initialisation of data sets
    for i in [ 1 .. a ] do
        A[i] := [];
        E[i] := [];
	K[i] := [];
        v[i] := [];
        w[i] := [];
        for k in [ 1 .. b ] do
            A[i][k] := [];
            E[i][k] := [];
            v[i][k] := [];
            w[i][k] := [];
        od;
        for h in [ 1 .. a ] do
            K[i][h] := [];
        od;
    od;
    for i in [ 1 .. b ] do
        M[i] := [];
        F[i] := [];
        for k in [ 1 .. a ] do
            M[i][k] := [];
            F[i][k] := [];
        od;
    od;
    for k in [ 1 .. b ] do
        D[k] := rec(remnant := [], pivots := []);
        B[k] := [];
        R[k] := [];
        X[k] := [];
        for j in [ 1 .. b ] do
            B[k][j] := [];
            R[k][j] := [];
            X[k][j] := [];
        od;
    od;
  
    # Step1: ClearDown
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            if j = 1 then       
               E[i][1] := 0*[ 1 .. DimensionsMat(Inp[i][j])[1] ];
            fi;

            tmp := ClearDown( f,C[i][j],D[j].pivots,D[j].remnant );
            A[i][j] := tmp[3];
            D[j].remnant := tmp[1]; 
            D[j].pivots := tmp[2];
           
            if not j = 1 then 
                E[i][j] := ExtendPivotRows( E[i][j-1],tmp[3][5] );
                v[i][j] := MKR( E[i][j-1],E[i][j] );
                w[i][j] := MKw( E[i][j-1],E[i][j] );
            else
                v[i][1] := E[i][1];
                E[i][1] := ExtendPivotRows( E[i][1],tmp[3][5] );
                v[i][1] := MKR(v[i][1],E[i][1] );
                w[i][j] := E[i][1];
            fi;

            for k in [ j+1 .. b ] do
                tmp := UpdateRow( f,A[i][j],C[i][k],B[j][k] );
                C[i][k] := tmp[1];
                B[j][k] := tmp[2];
            od;
        
            for h in [ 1 .. i-1 ] do 
                tmp := UpdateTrafo( f,A[i][j],K[i][h],M[j][h],v[h][j],0,w[h][j] );
                K[i][h] := tmp[1];
                M[j][h] := tmp[2];
            od;
            tmp := UpdateTrafo( f,A[i][j],K[i][i],M[j][i],v[i][j],1,w[i][j] );
            K[i][i] := tmp[1];
            M[j][i] := tmp[2];
            #Error( "DBUG----LOCATION1");
        od;
        
        # produce F
        for j in [ 1 .. b ] do
            F[j][i] := MKR( E[i][j],E[i][b] );
        od;
    od;
    #Error("DBUG----LOCATION2");
    
    # Step2: Riffle missing collumns into the M_jh's
    for j in [ 1 .. b ] do
        for h in [ 1 .. a ] do
            M[j][h] := RowLenghten( f,M[j][h],F[j][h] );
        od;
    od;

    # Step3: Upwards Cleaning
    for i in [ 1 .. b ] do #we use i for b-k+1 from now on
       k := b-i+1;
       R[k][k] := ImmutableMatrix(f,ShallowCopy(D[k].remnant));
    od; 
    for i in [ 1 .. b ] do
        k := b-i+1;
        for j in [ 1 .. k-1 ] do
           tmp := CEX( f,BitstringToCharFct(D[k].pivots),B[j][k] );
           X[j][k] := ImmutableMatrix(f,tmp[1]);
           R[j][k] := ImmutableMatrix(f,tmp[2]);
           #Error("check X"); 
           for h in [ k .. b ] do
              if not IsEmpty(R[k][h]) then
                 R[j][h] := R[j][h] + X[j][k]*R[k][h];
              fi;
           od;
           for h in [ 1 .. a ] do
              if not IsEmpty(M[k][h]) then
                 M[j][h] := M[j][h] + X[j][k]*M[k][h];
              fi;
           od;
        od;
    od;

    ## Write output

    # Beginning with row-select bitstring now named v
    v := [];
    rank := 0;
    for i in [ 1 .. a ] do
        v := Concatenation( v,E[i][b] );
        rank := rank + Sum( E[i][b] );
    od;

    # Glueing the blocks of M
    B := NullMat( rank,Length(v),f );
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
        for i in [ 1 .. a ] do
            if IsEmpty(M[j][i]) then
                M[j][i] := NullMat( rows[j],Length(E[i][b]),f );
            else
                M[j][i] := TransposedMat( RRF( NullMat( Length(E[i][b])-DimensionsMat(M[j][i])[2],DimensionsMat(M[j][i])[1],f ),TransposedMat(M[j][i]),E[i][b] ) );
            fi;

            B{[tmpR .. tmpR + DimensionsMat(M[j][i])[1]-1 ]}{[tmpC .. tmpC + DimensionsMat(M[j][i])[2]-1 ]}
            := M[j][i];
            tmpC := tmpC + DimensionsMat(M[j][i])[2];
        od;
        tmpR := tmpR + rows[j];
    od;

    # Glueing the blocks of R
    #Error("BREAKPOINT - Checking reamnant"); 
    C := NullMat( rank,ncols-rank,f );
    rows := [];
    w := [];

    for i in [ 1 .. b ] do
         rows[i] := 0;
         if IsEmpty(D[i].pivots) then
             w := Concatenation( w,0*[1..DimensionsMat(Inp[1][i])[2]] );
         else
             w := Concatenation( w,D[i].pivots );
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
                 if not IsEmpty(D[j].pivots) then
                     tmpC := tmpC + Sum( 1 - D[j].pivots );
                 elif  not IsEmpty(R[1][j]) then
                     tmpC := tmpC + DimensionsMat(R[1][j])[2];
                 fi;
                 continue;            
             fi;
 
             C{[tmpR .. tmpR + DimensionsMat(R[i][j])[1]-1 ]}{[tmpC .. tmpC + DimensionsMat(R[i][j])[2]-1 ]}
             := R[i][j];
             tmpC := tmpC + DimensionsMat(R[i][j])[2];
         od;
         tmpR := tmpR + rows[i];
     od;    

     # SLOW - only for testing 
     C := TransposedMat( RRF( TransposedMat(C), -IdentityMat( rank,f ),w  ) );

     return [v,w,C,B,D];
end;


