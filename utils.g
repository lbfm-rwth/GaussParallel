CEX := function( f,P,C )
 local i,L,R,col,dim,numr,numc,Lct,Rct;

 if C = [] then return [[],[]]; fi;

 dim := DimensionsMat( C );
 numr:=dim[1]; numc:=dim[2];
 Lct:=1;Rct:=1;
 L := NullMat(Size(P),numr,f); R:=NullMat(numc-Size(P),numr,f);
 for i in [1..numc] do
  col := TransposedMat(C)[i];
  if i in P then
   L[Lct] := col; Lct := Lct +1;
  else 
   R[Rct] := col; Rct := Rct +1;
  fi;
 od;
return [TransposedMat(L),TransposedMat(R)];
end;

BitstringToCharFct := function( P )
 local i,ct,Chi;
 ct:=1;
 Chi:=[];
 for i in [1..Length(P)] do
  if P[i] = 0 then ct := ct +1; continue; fi;
  ct := ct + 1; Add(Chi,i);
 od;
 return Chi;
end;

CharFctToBitstring := function( l,P )
 local i,Chi;
 Chi:=[];
 for i in [1..l] do
  if i in P then 
   Chi[i]:=1;
  else
   Chi[i]:=0;
  fi;
 od;
 return Chi;
end;

# TODO
ComplementBitstring := function( bstring )
end;

REX := function( f,T,C )
 local ret;
 
 if C = [] then return [[],[]]; fi;
 
 ret := CEX( f,T,TransposedMat(C) );
 return [TransposedMat(ret[1]),TransposedMat(ret[2])];
end;

PVC := function( s,t )
 local stOrdered,u,i;
 stOrdered := Concatenation(s,t); u:=0*[1..Length(s)+Length(t)];
 Sort(stOrdered);
 for i in t do
  u[Position(stOrdered,i)]:=1;
 od;
  
 return [stOrdered,u];
end;

RRF := function( R,RR,u )
 local ind,indR,indRR,Rnew;
 indR:=1; indRR:=1;
 Rnew := ShallowCopy(R);

 if R = [] then return RR; fi; if RR = [] then return R; fi;

 while (indR <= DimensionsMat(R)[1]) or (indRR <= DimensionsMat(RR)[1]) do
  ind := indR + indRR -1;
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

MCP := function( X )
 return ShallowCopy( X );
end;

ECH := function( f,H )
 local EMT,m,k,M,K,R,S,N,r,s,t,i,ind,Id,one,zero;

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
 
 ind := 1;
 for i in [1..DimensionsMat(m)[1]] do
  if m[i] = zero*[1..DimensionsMat(m)[2]] then 
   Add(s,0);  
  else
   Add(M,m[i]);
   if not k = [] then
    Add(K,k[i]);
   fi;
   Add(s,1);
  fi;
 od;
 
 ind := 1;
 for i in [1..DimensionsMat(r)[1]] do
  if ind > DimensionsMat(Id)[1] then
   Add(t,0);
   Add(R,r[i]);
   continue;
  fi;
  if r[i] = Id[ind] then
   Add(t,1);
   ind := ind + 1;
  else
   Add(t,0);
   Add(R,r[i]);
  fi;
 od; 
 
 return [-TransposedMat(M),TransposedMat(K),-TransposedMat(R),s,t];
end;

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


