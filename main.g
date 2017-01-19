### TODO ###
## rename variables, e.g.:
## t -> ListOfPivotColumns
### TODO ###
###############################
# function ClearColumn
# Input:
#   f - Field
#   H - Input matrix
#   t - Bitlist of pivot columns positions
#   R - Residue of above block after echolonization
#
# Output:
#   [ R, t, T ]
#   R - ??
#   t - ??
#   T - ??
###############################
ClearColumn := function( f, H, t, R )
    local tmp, Chi, ct, HH, tt, RR, ttt, RRR, i, RRn, A, AA, T, M, K, E, s, u;
    if not Length(t) = DimensionsMat(H)[2] then
        Error( "Length of bitlist t does not match dimensions of H!" );
    fi;

    #### INITIALIZATION ####
    ## Residue R was empty. Thus matrix above has full column-rank.
    if IsEmpty( R ) then
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

    ## First step (i=1,  i block row) or all matrices above have rank 0
    if IsEmpty( t ) then
        # A is empty iff this case happens
        A := [];
        HH := H;
    else
        # Column Extraction
        tmp := CEX( f, BitstringToCharFct(t), H );
        A := tmp[1]; AA := tmp[2];
        # Reduce H to (0|HH)
        # Mult Add
        HH := AA + A*R;
    fi;
    #### END INITIALIZATION ####

    # Echelonization
    tmp := ECH( f, HH );
    M:=tmp[1];K:=tmp[2];RR:=tmp[3];s:=tmp[4];tt:=tmp[5];
    Error( "Break Point - echel" );

    # TODO complement then extend?
    Chi := [];ct:=1;
    if not tt=[] then
        for i in [1..Length(t)] do
            if t[i]=0 then
                if tt[ct]=1 then Add(Chi, i); fi;
                ct:= ct+1;
            fi;
        od;
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
        return [RR, ttt, T];
    fi;
    # If RR is empty, but tt is not, then the bitstring tt, representing
    # the positions of the new pivot columns, is AllOne.
    # In this case, there is nothing to be done here.

    Error( "Break Point - before CEX new residue" );
    tmp := CEX( f, BitstringToCharFct(tt), R );
    E:=tmp[1];RRn:=tmp[2];
    ## Update the residue and the pivot column bitstring
    tmp := PVC( BitstringToCharFct(t), Chi );
    ttt:=CharFctToBitstring(Length(t), tmp[1]); u:=tmp[2];
    Error( "Break Point - after CEX new residue" );

    ## Did column extraction return empty values?
    ## E should never be empty in this context!
    ## TODO proof?
    if IsEmpty(E) then
        Error( "This should not have happened!" );
    fi;
    ## RRn is empty, iff. the new pivot columns completely
    ## annihilate the old residue.
    if IsEmpty(RRn) then
        RR := [];
    else
        RRR:=RRn+E*RR;
        RR := RRF( RRR, RR, u );
    fi;

    T:=[A, M, E, K, s, u];
    return [RR, ttt, T];
end;

UpdateRow := function( f, T, H, Bjk )
 local A, E, M, K, s, u,  tmp, Z, V, X, W, S, B;
 B := Bjk;
 A:=T[1];M:=T[2];E:=T[3];K:=T[4];s:=T[5];u:=T[6];
 
 ###
 # If A is empty, there are no rowoperations form above to consider
 ###
 if IsEmty(A) then
  Z := H:
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
  H :=K*V+W;
 fi;

 return [H, B];
end;

BackClean := function( k, Bjk, T, Cij )
    return 0;
end;

Step1 := function( A )
 local C, n, f, tmp,   i, j, k, B, T, Rj, H, tj, V, W;
 f := DefaultFieldOfMatrix( A );
 n := 2; #Gcd( DimensionsMat(A)[1], DimensionsMat(A)[2] );
 C := ChopMatrix( n, A );
 B := ShallowCopy(C); #Init B as nxn
 Rj := [];
 tj := [];

 for i in [1..n] do
  for j in [1..n] do
   H := C[i][j];

   if i = 1 then
    tmp := ECH( f, H );
    tj[j] := tmp[5];
    Rj[j] := tmp[3];
    T := [
        [],
        tmp[1],
        CEX( f, BitstringToCharFct(tj[j]), Rj[j] )[1],
        tmp[2],
        tmp[4],
        tj[j]
    ];
   else
    tmp := ClearColumn( f, H, tj[j], Rj[j] );
    Rj[j] := tmp[1];
    tj[j] := tmp[2];
    T := tmp[3];
   fi;

   for k in [j+1..n] do
    if i = 1 then
     tmp := REX( f, BitstringToCharFct(T[5]), C[i][k] );
     V:=tmp[1];
     if V = [] then
      B[j][k]:=[];
      continue;
     fi;
     B[j][k] := ShallowCopy(T[2]*V);
    fi;

    tmp := UpdateRow( i, f, T, C[i][k], B[j][k] );
    C[i][k] := tmp[1];
    B[j][k] := tmp[2];
   od;

  od;
 od;

 return [C, B, Rj, tj];
end;
