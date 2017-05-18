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
        HH := ImmutableMatrix(f,AA) + ImmutableMatrix(f,A)*R;
    fi;
 
    # Echelonization
    tmp := ECH( f, HH );
    M:=ImmutableMatrix(f,tmp[1]);K:=ImmutableMatrix(f,tmp[2]);RR:=ImmutableMatrix(f,tmp[3]);s:=tmp[4];tt:=tmp[5];

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
        #  RR := R; # Residue does not change
        #  ttt := t;
        u := 0 * t;
        T:=[ImmutableMatrix(f,A), M, E, K, s, u];
        return [R, t, T ];
    fi;
    # If RR is empty, but tt is not, then the bitstring tt, representing
    # the positions of the new pivot columns, is 1 everywhere.
    # In this case, there is nothing to be done here.

    tmp := CEX( f, BitstringToCharFct(tt), R );
    E:=ImmutableMatrix(f,tmp[1]); RRn:=ImmutableMatrix(f,tmp[2]);
    ## Update the residue and the pivot column bitstring
    tmp := PVC( BitstringToCharFct(t), BitstringToCharFct(Chi) );
    ttt:=CharFctToBitstring(DimensionsMat(H)[2], tmp[1]); u:=tmp[2];
    
    T:=[ImmutableMatrix(f,A), M, E, K, s, u];
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
 V:=ImmutableMatrix(f,tmp[1]);W:=ImmutableMatrix(f,tmp[2]);
 ###
 # If V is empty, then there where no operations exept from A
 # in this case there is nothing more to update
 ###
 if IsEmpty(V) then
  return [Z,B]; 
 else 
  X:=M*V;
 fi;

 if IsEmpty(E) then
  S:= [];
 else
  S:= E*X+B;
 fi;
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

  return ImmutableMatrix(f,TransposedMat( new ));
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
      V := ImmutableMatrix(f,tmp[1]);
      W := ImmutableMatrix(f,tmp[2]);
   else
      if not IsEmpty(A) then
         Z := Y + A*M;
      else
         Z := Y + M;
      fi;         
      
      tmp := REX( f,BitstringToCharFct(s),Z );
      V := ImmutableMatrix(f,tmp[1]);
      W := ImmutableMatrix(f,tmp[2]);
   fi;

   if not IsEmpty(V) then
       X := MM*V;
       if IsEmpty(E) then
        S:= [];
       else
        S := M + E*X;     
       fi;
       M := RRF( S,X,u );
   fi;
   
   #Error("ende");      
   if not IsEmpty( V ) and not IsEmpty( KK ) then
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

ChoppedMatrix := function( f,A,nrows,ncols )
    local i,j,rrem,crem,AA,a,b;
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

RowLenghten := function( f,M,F )
   local new,i,current;
   if IsEmpty(M) then
      return M;
   fi;
 
   current := 1; 
   new := NullMat( Length(F),DimensionsMat(M)[1],f );
   for i in [ 1 .. Length(F) ] do
       if F[i] = 1 then
          new[i] := TransposedMat(M)[current]; 
          current := current + 1;
       fi;
   od;
   return ImmutableMatrix(f,TransposedMat(new));
end;
