#####################################################################
# INPUT:  f - a finite field,
#         H - a Matrix over f,
#         t - a bitstring containing the pivotal
#             columns of the current block-column,
#         R - a Matrix over f, containing the non-pivotal columns.
# OUTPUT: the new remnant together with the updated bitstring 
#         and the exput (see paper). 
#####################################################################
ClearDown := function( f, H, t, R )
    local tmp, Chi, ct, HH, QQ, Q1,Q2,Q3, tt, RR, ttt, RRR, i, RRn, A, AA, T, M, K, E, s, u;
   
    if IsEmpty(H) then
        return [R,t,[[],[],[],[],[],[]]];
    fi; 

    ## Residue R was empty. Thus matrix above has full column-rank.
    if not IsEmpty(t) and IsEmpty( R ) then
        A := H;
        M := [];
        E := [];
        K := [];
        s := [];
        RR := [];
        tt := [];
        ttt := t;
        u := 0 * t;
        T:=[A, M, E, K, s, u];
        return [RR, ttt, T];
    fi;

    ## First step (first block row) or all matrices above have rank 0
    if IsEmpty( t ) then
        # A is empty iff this case happens
        A := [];
        HH := H;
    else
        # Column Extraction
        tmp := CEX( f, BitstringToCharFct(t), H );
        A := tmp[1]; AA := tmp[2];
        # Reduce H to (0|HH)
        HH := AA + A*R;
    fi;
 
    # Echelonization
    tmp := ECH( f, HH );
    M:=tmp[1];K:=tmp[2];RR:=tmp[3];s:=tmp[4];tt:=tmp[5];

    ## Producing a riffle -Chi- from old and new pivot columns
    if IsEmpty(t) then Chi := tt;
    else
     Chi := 0*[1..DimensionsMat(H)[2]];ct:=1;
     if not tt=[] then
        for i in [1..Length(Chi)] do
            if t[i]=0 then
                if tt[ct]=1 then Chi[i] := 1; fi;
                ct:= ct+1;
            fi;
        od;
     fi;
    fi;

    ## Are we running into any special cases, where return values
    ## of the above echelonization are empty?
    # The case K empty is handled by UpdateRow.
    #
    # The following are equivalent:
    # s is empty
    # tt is empty
    # M is empty
    #
    # This only happens when A*R + AA = 0
    if ForAny( [ M, s, tt ], IsEmpty ) then
        M := [];
        E := [];
        K := [];
        s := [];
        RR := R; # Residue does not change
        ttt := t;
        u := 0 * t;
        T:=[A, M, E, K, s, u];
        return [RR, ttt, T ];
    fi;
    # If RR is empty, but tt is not, then the bitstring tt, representing
    # the positions of the new pivot columns, is 1 everywhere.
    # In this case, there is nothing to be done here.

    tmp := CEX( f, BitstringToCharFct(tt), R );
    E:=tmp[1];RRn:=tmp[2];
    ## Update the residue and the pivot column bitstring
    tmp := PVC( BitstringToCharFct(t), BitstringToCharFct(Chi) );
    ttt:=CharFctToBitstring(DimensionsMat(H)[2], tmp[1]); u:=tmp[2];
    
    T:=[A, M, E, K, s, u];
    if IsEmpty(E) then
        return [RR,ttt,T];
    fi;

    ## Did column extraction return empty values?
    ## RRn is empty, iff. the new pivot columns completely
    ## annihilate the old residue.
    if IsEmpty(RRn) then
        RR := [];
    else
        RRR:=RRn+E*RR;
        RR := RRF( RRR, RR, u );
    fi;
     
    return [RR, ttt, T ];
end;

#####################################################################
# INPUT:  f - a finite field,
#         H - a matrix tracking the changes on the pivot rows 
#             of the input-matrix
#         T - the exput from previous ClearDowns,
#         Bjk - a matrix tracking the changes on the non-pivot rows of 
#               the input-matrix.
# OUTPUT: the updated matrices H and Bjk according to new found pivot rows
#         tracked by T.
#####################################################################
UpdateRow := function( f, T, H, Bjk )
 local A, E, M, K, s, u,  tmp, Z, V, X, W, S, B;
 B := Bjk;
 A:=T[1];M:=T[2];E:=T[3];K:=T[4];s:=T[5];u:=T[6];
 
 ###
 # If A is empty, there are no rowoperations form above to consider
 ###
 if IsEmpty(A) then
  Z := H;
 else 
  Z := A*B+H;
 fi;

 tmp := REX( f, BitstringToCharFct(s), Z );
 V:=tmp[1];W:=tmp[2];
 ###
 # If V is empty, then there where no operations exept from A
 # in this case there is nothing more to update
 ###
 if IsEmpty(V) then
  return [Z,B]; 
 else 
  X:=M*V;
 fi;

 S:= E*X+B;
 B:=RRF( S, X, u );
 
 ###
 # if K is empty, then s is the all-one-bitstring and there are no non-pivot rows
 # which would change according to K. So K should be empty and there is nothing more to update
 ###
 if not IsEmpty(K) then
  # s is neither empty nor all-one at this point
  H := K*V+W;
 else
  H := W;
 fi;
 

 return [H, B];
end;

RiffleIn := function( X,u,w,type,f )
  local i,F,new;
  
  if 1-u = 0*[1..Length(u)] then
      return X;
  fi;

  if IsEmpty(X) then
      F := IdentityMat( Length(w),f );
  else    
      F := IdentityMat( DimensionsMat(X)[1],f );
  fi;
  if type = 0 then
      F := Zero(f) * F;
  fi;
  F := TransposedMat( CEX( f,BitstringToCharFct(w),F )[1] );  
  new := RRF( F,TransposedMat(X),u );

  return TransposedMat( new );
end;

#####################################################################
# INPUT:  f - a finite field,
#         M - a matrix tracking the changes on the pivot rows 
#             of the input-matrix
#         T - the exput from previous ClearDowns,
#         K - a matrix tracking the changes on the non-pivot rows of 
#               the input-matrix.
#         v,flag,w -?
# OUTPUT: the updated matrices M and K according to new found pivot rows
#         tracked by T.
#####################################################################
UpdateTrafo := function( f, T, K, M,v,flag,w )
 local A, E, MM, KK, s, u,  i,tmp, riffle, Z, V, X, W, Y,S;
 A:=T[1];MM:=T[2];E:=T[3];KK:=T[4];s:=T[5];u:=T[6];
   #Error("anfang"); 
   Y := RiffleIn( K,v,w,flag,f ); 
 
   if IsEmpty(A) and IsEmpty(MM) then
       return [K,M];
   fi;

   if IsEmpty(M) then

     #if  IsEmpty(K) then
     #       if  IsEmpty(MM) then
     #           return [K,M];
     #       fi;
     #       X := MM;
     #       S := E*MM;
     #       M := RRF( S,X,u );
     #       #Error( "richtig?" );
     #       return [KK,M];
     # fi; 

      tmp := REX( f,BitstringToCharFct(s),Y );
      V := tmp[1];
      W := tmp[2];
   else
      if not IsEmpty(A) then
         Z := Y + A*M;
      else
         Z := Y + M;
      fi;         
      
      tmp := REX( f,BitstringToCharFct(s),Z );
      V := tmp[1];
      W := tmp[2];
   fi;

   if not IsEmpty(V) then
       X := MM*V;
       S := M + E*X;     
       M := RRF( S,X,u );
   fi;
   
   #Error("ende");      
   if not IsEmpty( V ) then
      K := W + KK*V;
   else
      K := W;
   fi;

   return [K, M];
end;

ExtendPivotRows := function( old,new )
    local current,i; current := 1;
    old := MutableCopyMat(old);
    for  i in [ 1 .. Length(old) ] do
        if current>Length(new) then
            return old;
        fi;
        if old[i]=0 then
            if new[current]=1 then
                old[i] := 1;
            fi;
            current := current +1;
        fi;
    od;
    return old;
end;

ChoppedMatrix := function( A,nrows,ncols )
    local i,j,rrem,crem,AA,a,b;
    rrem := DimensionsMat(A)[1] mod nrows;
    crem := DimensionsMat(A)[2] mod ncols;
    a := ( DimensionsMat(A)[1] - rrem ) / nrows; 
    b := ( DimensionsMat(A)[2] - crem ) / ncols; 
    AA := [];

    for  i  in [ 1 .. nrows-1] do
        AA[i] := [];
        for j in [ 1 .. ncols-1 ] do
            AA[i][j] := A{[(i-1)*a+1 .. i*a]}{[(j-1)*b+1 .. j*b]};
        od;
    od;
    AA[nrows] := [];
    for i in [ 1 .. nrows-1 ] do
        AA[i][ncols] := A{[(i-1)*a+1 .. i*a]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]};
    od;
    for j in [ 1 .. ncols-1 ] do
        AA[nrows][j] := A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(j-1)*b+1 .. j*b]};
    od;
    AA[nrows][ncols] := A{[(nrows-1)*a+1 .. DimensionsMat(A)[1]]}{[(ncols-1)*b+1 .. DimensionsMat(A)[2]]};    return AA;
end;

MKR := function( interm,final )
    local new,i;
    new := [];
 
    for i in [ 1 .. Length(final) ] do
        if final[i] = 0 then continue; fi;
        if interm[i] = 1 then 
     	    Add( new,1 );
        else
            Add( new,0 );
        fi;
    od;
    return new;
end;

MKw := function( old,new )
    local tmp,i;
    tmp := [];
 
    for i in [ 1 .. Length(old) ] do
        if new[i] = 0 then
            Add( tmp,0 );
        fi;
        if new[i] = 1 and old[i] = 0 then
            Add( tmp,1 );
        fi;
    od;
    return tmp;
end;

RowLenghten := function( M,F )
   local new,i,current;
   if IsEmpty(M) then
      return M;
   fi;
 
   current := 1; 
   new := NullMat( Length(F),DimensionsMat(M)[1],DefaultFieldOfMatrix(M) );
   for i in [ 1 .. Length(F) ] do
       if F[i] = 1 then
          new[i] := TransposedMat(M)[current]; 
          current := current + 1;
       fi;
   od;
   return TransposedMat(new);
end;

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
            M[j][h] := RowLenghten( M[j][h],F[j][h] );
        od;
    od;

    # Step3: Upwards Cleaning
    for i in [ 1 .. b ] do #we use i for b-k+1 from now on
       k := b-i+1;
       R[k][k] := ShallowCopy(D[k].remnant);
    od; 
    for i in [ 1 .. b ] do
        k := b-i+1;
        for j in [ 1 .. k-1 ] do
           tmp := CEX( f,BitstringToCharFct(D[k].pivots),B[j][k] );
           X[j][k] := tmp[1];
           R[j][k] := tmp[2];
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


