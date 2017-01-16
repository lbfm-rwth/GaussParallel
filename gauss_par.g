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

ClearColumn := function( f,H,t,R )
 local tmp,Chi,ct,HH,tt,RR,ttt,RRR,i,RRn,A,AA,T,M,K,E,s,u;

 if H = [] then return [ [],[],[ [],[],[],[],[] ] ]; fi;
 
 if R = [] then 
  u := []; for i in [1..DimensionsMat(H)[1]] do Add(u,1); od;
  return [R,t, [H,[],[],[],[],u] ]; 
 fi; 

 tmp := CEX( f,BitstringToCharFct(t),H );
 A := tmp[1]; AA := tmp[2];
 HH := AA + A*R;
 tmp := ECH( f,HH ); 
 M:=tmp[1];K:=tmp[2];RR:=tmp[3];s:=tmp[4];tt:=tmp[5];
 
 Chi := [];ct:=1;
 if not tt=[] then
  for i in [1..Length(t)] do
   if t[i]=0 then  
    if tt[ct]=1 then Add(Chi,i); fi;
    ct:= ct+1;
   fi;   
  od;
 fi;

 tmp := CEX( f,BitstringToCharFct(tt),R );
 E:=tmp[1];RRn:=tmp[2];
  
 if E = [] then
  u := []; for i in [1..DimensionsMat(H)[1]] do Add(u,1); od; 
  return [R,t, [A,M,E,K,s,u] ];
 fi;
 if RRn = [] then RRR := Zero(f)*RR; else  
  RRR:=RRn+E*RR;
 fi;

 tmp := PVC( BitstringToCharFct(t),Chi );
 ttt:=CharFctToBitstring(Length(t),tmp[1]);u:=tmp[2];
 RR := RRF( RRR,RR,u ); 
 
 T:=[A,M,E,K,s,u];
 return [RR,ttt,T];
end;

UpdateRow := function( i,f,T,H,Bjk )
 local A,E,M,K,s,u, tmp,Z,V,X,W,S,B;
 B := Bjk;
 A:=T[1];M:=T[2];E:=T[3];K:=T[4];s:=T[5];u:=T[6];
 
 Z := A*B+H;
 tmp := REX( f,BitstringToCharFct(s),Z );
 V:=tmp[1];W:=tmp[2];
 X:=M*V;
 
 if i = 1 then return [K*V+W,X ]; fi;
 
 S:= E*X+B;
 B:=RRF( S,X,u ); 
 H :=K*V+W;
 
 return [H,B];
end;

BackClean := function( k,Bjk,T,Cij )
 return 0;
end;

Step1 := function( A )
 local C,n,f,tmp,  i,j,k,B,T,Rj,H,tj,V,W;
 f := DefaultFieldOfMatrix( A );
 n := 2; #Gcd( DimensionsMat(A)[1],DimensionsMat(A)[2] );
 C := ChopMatrix( n,A );
 B := ShallowCopy(C); #Init B as nxn
 Rj := [];
 tj := [];

 for i in [1..n] do
  for j in [1..n] do
   H := C[i][j];
   
   if i = 1 then
    tmp := ECH( f,H );
    tj[j] := tmp[5];
    Rj[j] := tmp[3];  
    T := [ [],tmp[1], CEX( f,BitstringToCharFct(tj[j]),Rj[j] )[1], tmp[2],tmp[4],tj[j] ];   
   else
    tmp := ClearColumn( f,H,tj[j],Rj[j] );
    Rj[j] := tmp[1];
    tj[j] := tmp[2];
    T := tmp[3];
   fi;

   for k in [j+1..n] do
    if i = 1 then
     tmp := REX( f,BitstringToCharFct(T[5]),C[i][k] );
     V:=tmp[1];
     if V = [] then 
      B[j][k]:=[];
      continue; 
     fi;
     B[j][k] := ShallowCopy(T[2]*V);
    fi;
     
    tmp := UpdateRow( i,f,T,C[i][k],B[j][k] );  
    C[i][k] := tmp[1];
    B[j][k] := tmp[2];
   od;   

  od;
 od;

 return [C,B,Rj,tj];
end;



















##
#It follows a sequencial version of the parallel gauss alg
##
gaussElem := function( XX )
 local f,chsz,a,b,tmp,i,j,k,m,C,  Aij,D,Eij,G,H,Kij,L,Mij,N,Q,Uij,S,Tij,V,W,X,Z,P,R,B;
 chsz := Gcd( DimensionsMat(XX)[1],DimensionsMat(XX)[2] );
 b := DimensionsMat(XX)[1]/chsz;
 a := DimensionsMat(XX)[2]/chsz;

 P:=[]; R:=[]; B:=[];
 f := DefaultFieldOfMatrix( XX );
 C := ChopMatrix( chsz,XX );
 
 for i in [1..chsz] do
  B[i]:=[]; 
  for j in [1..chsz] do
  B[i][j] :=[]; 
  od;
 od;

 for i in [1..chsz] do
 
  for j in [1..chsz] do
 
   if i = 2 then return [C,B,R,P]; fi;
   
   if i = 1 then
    H := ShallowCopy(C[i][j]);
    tmp := ECH( f,H );
    Mij:=tmp[1];Kij:=tmp[2];R[j]:=tmp[3];Tij:=tmp[4];P[j]:=tmp[5];
   else
    tmp := CEX( f,BitstringToCharFct(P[j]),C[i][j] );
    Aij:=tmp[1];G:=tmp[2];
  
    if not Aij=[] and not R[j]=[] and not G=[] then
     H := Aij*R[j]+G;
    else H:=[]; fi;
  
    tmp := ECH( f,H );
    Mij:=tmp[1];Kij:=tmp[2];N:=tmp[3];Tij:=tmp[4];Q:=tmp[5];

    tmp := CEX( f,BitstringToCharFct(Q),R[j] );
    Eij:=tmp[1];D:=tmp[2];
  
    if not Eij=[] and not N=[] and not D=[] then
    L := Eij*N+D; else L:=[]; fi;
 
    tmp := PVC( BitstringToCharFct(P[j]),BitstringToCharFct(Q) );
    P[j]:=CharFctToBitstring(a,tmp[1]);Uij:=tmp[2];
 
    tmp := RRF( N,L,Uij ); #or N,L,Uij ??
    R[j] := tmp;
   fi;

   for k in [j+1..chsz] do
    m := j+1;
    if i = 1 then
     tmp := REX( f,BitstringToCharFct(Tij),C[i][k] );
     V:=tmp[1];W:=tmp[2];
     
     if not Mij=[] and not V=[] then
     X := Mij*V; else X:=[]; fi;
  
     B[j][k] := ShallowCopy(X);
    else
     if not Aij=[] and not B[j][k]=[] and not C[i][k]=[] then
     Z := Aij*B[j][k]+C[i][k]; else Z:=[]; fi;
    
     tmp := REX( f,BitstringToCharFct(Tij),Z );
     V:=tmp[1];W:=tmp[2];
     
     if not Mij=[] and not V=[] then
     X := Mij*V; else X:=[]; fi;
  
     if not Eij=[] and not X=[] and not B[j][k]=[] then
     S := Eij*X+B[j][k]; else S:=[]; fi;
  
     B[j][k] := RRF( S,X,Uij );
    fi;
   
    if m<chsz then
     if not Kij = [] and not V=[] and not W=[] then
     C[i][k]:=Kij*V+W; fi;
    fi;
   od; 
  od;
 od; 

 return [C,B,R,P];
end;











