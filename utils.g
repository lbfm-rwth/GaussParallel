CEX := function( f,P,C )
 local i,L,R,numc,col,Lct,Rct;

 if IsEmpty(C) then return [[],[]]; fi;

 numc := DimensionsMat( C )[2];
 C := TransposedMat(C);
 Lct:=1; Rct:=1;
 L := []; R := []; 
 for i in [1..numc] do
  col := C[i];
  if i in P then
   L[Lct] := col; Lct := Lct +1;
  else 
   R[Rct] := col; Rct := Rct +1;
  fi;
 od;
return [TransposedMat(L),TransposedMat(R)];
end;

BitstringToCharFct := function( P )
 local i,Chi,l,insertCT;
 insertCT := 1;
 Chi:=[];
 l := Length(P);
 for i in [1..l] do
  if P[i] = 0 then continue; fi; 
  Chi[insertCT] := i;
  insertCT := insertCT + 1;
 od;
 return Chi;
end;

CharFctToBitstring := function( l,P )
 local i,Chi;
 Chi := [];
 for i in [ 1 .. l ] do
     Chi[i] := 0;
 od;
 for i in P do 
   Chi[i]:=1;
 od;
 return Chi;
end;

REX := function( f,TT,C )
 local ret;
 
 if IsEmpty(C) then return [[],[]]; fi;
 ret := CEX( f,TT,TransposedMat(C) );
 return [TransposedMat(ret[1]),TransposedMat(ret[2])];
end;

PVC := function( s,t )
 local stOrdered,u,i,l;
 stOrdered := Concatenation(s,t); 
 l := Length(s)+Length(t);
 u := [];
 for i in [ 1 .. l ] do
     u[i] := 0;
 od;
 Sort(stOrdered);
 for i in t do
  u[Position(stOrdered,i)]:=1;
 od;
  
 return [stOrdered,u];
end;

RRF := function( R,RR,u )
 local l,ind,indR,indRR,Rnew;
 indR:=1; indRR:=1;
 Rnew := [];
 l := Length(u);
 if R = [] then return RR; fi; if RR = [] then return R; fi;

 while (indR <= DimensionsMat(R)[1]) or (indRR <= DimensionsMat(RR)[1]) do
  ind := indR + indRR -1;
  if l < ind then
      break;
  fi;
  if u[ind] = 0 then
   Rnew[ind] := R[indR];
   indR := indR + 1;
  else
   Rnew[ind] := RR[indRR];
   indRR := indRR + 1;
  fi;
 od;
  
 return Rnew;
end;

ECH := function( f,H )
 local sct,Mct,Kct,tct,Rct,EMT,m,k,M,K,R,S,N,r,s,t,i,ind,Id,one,zero,dims,dimId;

 if H = [] then return [ [],[],[],[],[] ]; fi;
 if Rank(H) = 0 then return [ [],[],[],[],[] ]; fi;
 
 EMT := EchelonMatTransformation( H );
 m := TransposedMat(EMT.coeffs);
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
   if not k = [] then
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
 
 return [-TransposedMat(M),TransposedMat(K),-TransposedMat(R),s,t];
end;


#####
# The following functions are not in use in current versions - use careful
#####
testECH := function( A )
 local f,indS,i,L,R,indSS ,Id,S,N,SS,NN;
 f := DefaultFieldOfMatrix( A );
 L := ECH( f,A );
 R := EchelonMatTransformation( A );
 
 Id := IdentityMat( Length(R.vectors),f );
 indS := 1; indSS := 1;
 S := []; SS := []; N:=[]; NN:=[];
 #calculate S,N,S'N' from L
 for i in [1..Length(L[4])] do
  if L[4][i] = 1 then
   Add( S,Id[indS] ); 
   indS := indS + 1;
  else
   Add( S,(0*Id[1][1])*Id[1] );
  fi;
 od;
 for i in [1..Length(L[5])] do
  if L[5][i] = 1 then
   Add( SS,Id[indSS] ); 
   indSS := indSS + 1;
  else
   Add( SS,(0*Id[1][1])*Id[1] );
  fi;
 od;
 
 indS := 1;
 indSS := 1;
 if not DimensionsMat(S)[2] = DimensionsMat(A)[1] then
  Id := IdentityMat( DimensionsMat(A)[1] - DimensionsMat(S)[2],f );
  for i in [1..Length(L[4])] do
   if L[4][i] = 0 then
    Add( N,Id[indS] ); 
    indS := indS + 1;
   else
    Add( N,(0*Id[1][1])*Id[1] );
   fi;
  od;
 fi;
 if not DimensionsMat(SS)[2] = DimensionsMat(A)[2] then
  Id := IdentityMat( DimensionsMat(A)[2] - DimensionsMat(SS)[2],f );
  for i in [1..Length(L[5])] do
   if L[5][i] = 0 then
    Add( NN,Id[indSS] ); 
    indSS := indSS + 1;
   else
    Add( NN,(0*Id[1][1])*Id[1] );
   fi;
  od;
 fi;
 
 S := TransposedMat(S);
 SS := (SS);
 N := TransposedMat(N);
 NN := (NN);
 return [L[1]*S*A*SS,L[1]*S*A*NN,(L[2]*S+N)*A*SS,(L[2]*S+N)*A*NN,S,N,SS,NN];
end;

ChopMatrix := function( n,C )
 local b,a,i,j,k,RowChop,CChop;
 b := DimensionsMat( C )[1]/n;
 a := DimensionsMat( C )[2]/n;
 CChop := [[]];
 RowChop := [[]];
 
 for i in [1..n] do
  RowChop[i] := [];
  for j in [1..b] do 
   Add(RowChop[i],C[(i-1)*b+j]);
  od;
 od;
 
 for i in [1..n] do
  CChop[i]:=[];
  for k in [1..n] do
   CChop[i][k] := [];
   for j in [1..a] do
    Add(CChop[i][k],TransposedMat(RowChop[i])[(k-1)*a+j]);
   od;
    CChop[i][k] := TransposedMat(CChop[i][k]);
  od;
 od;
 
 return CChop;
end;

BuildPermutationMat := function( f,t )
 local Id,P,i,ct1,ct0;

 Id := IdentityMat( Length(t),f );
 P := [];
 ct1 := 1; ct0 := 1;
 for i in [1..Length(t)] do
  if t[i] = 1 then ct0 := ct0+1; fi;
 od; 
 for i in [1..Length(t)] do
  if t[i] = 1 then 
   P[i] := Id[ct1]; ct1 := ct1+1;
  else
   P[i] := Id[ct0]; ct0 := ct0+1;
  fi;
 od;
  
 return TransposedMat( P ); 
end;


