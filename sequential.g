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
