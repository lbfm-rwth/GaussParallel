### TODO ###
## rename variables, e.g.:
## t -> ListOfPivotColumns
### TODO ###

BackClean := function( f,k,n,B,t,C )
 local tmp,i,j,X,m;
 
 for i in [1..n] do
  for j in [1..n] do
   C[i][j] := MutableCopyMat(C[i][j]);
  od;
 od;
 for j in [1..k-1] do
  tmp := CEX( f,BitstringToCharFct( t ),B[j][k] );
  X := tmp[1]; C[j][k] := tmp[2];
  for m in [k..n] do
   if not IsEmpty( C[k][m] ) then
    C[j][m] := C[j][m] + X*C[k][m];
   fi;
  od;
 od; 
 
 return C;
end;

InsertBlock := function( f,X,t,R,C,i,j,ct,n,rct )
 local k,l,newct,B,Rct;
 newct:=0; Rct:=1; B := MutableCopyMat(X);
 if i=j then
  #if IsEmpty(R[j]) then return [B,ct]; fi;
 
  for k in [1..Length(t[j])] do
   if t[j][k] = 1 then
    B[ct+newct][Length(t[j])*(j-1)+k] := -One(f);
    newct := newct+1;
   else
    for l in [1..DimensionsMat(R[j])[1]] do
     B[ct-1+l][Length(t[j])*(j-1)+k] := R[j][l][Rct];
    od;
    Rct := Rct+1;
   fi;
  od;
 fi;
 
 if not i=j then
  if IsEmpty(C[i][j]) then return [B,ct]; fi;
  for k in [1..DimensionsMat(C[i][j])[1] ] do
   newct := 0;
   for l in [1..DimensionsMat(C[i][j])[2] ] do
    if not IsEmpty(t[j]) then
     while t[j][l+newct] = 1 do newct := newct+1; od;
    fi;
    B[rct-1+k][ (DimensionsMat(X)[2]/n)*(j-1)+l+newct ] := C[i][j][k][l];
   od;
  od;
  newct := 0;
 fi;

 return [B,ct+newct];
end;

Step2 := function( f,n,nrows,ncols,returnList )
 local M,Mj,MM,KK,shift,K,list,rank,i,j,S,s,current,Id,rct,ct,C,B,t,R,tmp, k,T;
 C := returnList[1];
 B := returnList[2];
 t := returnList[4];
 R := returnList[3];
 T := [];

 for k in [1..n] do
  if not IsEmpty(t[k]) then
   T := Concatenation( T,t[k] ); 
  else 
   T := Concatenation( T,0*[1..ncols/n] );
  fi;
  C[k][k] := R[k];
 od;

 for k in [1..n] do
  C := BackClean( f,n-k+1,n,B,t[n-k+1],C );
 od;
 
 #Error( "After BackClean" );
 ###
 # We now rearrange the results so that we can return the gauss normal form
 ### 
 rank := 0;
 for k in [1..Length(T)] do
  if T[k] = 1 then rank := rank +1; fi;
 od;
 
 B := MutableCopyMat(NullMat( rank,ncols,f ));
 Id := IdentityMat( rank, f );

 ## Build echelonform
 ct := 1;
 for k in [1..n] do
  rct := ct;
  for i in [k..n] do
   tmp := InsertBlock( f,B,t,R,C,k,i,ct,n,rct );
   B := tmp[1]; ct := tmp[2];
  od;
 od;

 return [B,BuildPermutationMat( f,T )];
end;

GaussParallel := function( A )
 local i,ct,k,f,n,l,rank,nrows,ncols;
 f := DefaultFieldOfMatrix( A );
 nrows := DimensionsMat(A)[1];
 ncols := DimensionsMat(A)[2];
 n := Gcd( nrows, ncols );
 
# Error( "Break Point - before Step1" );
 l := Step2( f,n,nrows,ncols, Step1( A,n ));
 return rec(vectors := -l[1], columnPermutation := l[2] );
end;

TestGaussParallel := function( nr,nc,iter )
 local i,test,A,bools,gap;
 
 bools:=[]; 
 for i in [1..iter] do
  A := RandomMat(3*nr,nc,GF(5));
  A[1] := Zero(GF(5))*A[1];
  A := RandomMat(nc,2*nr,GF(5)) * A;
  
  test := GaussParallel( A );
  gap := EchelonMat(A);

  bools[i] := test.vectors = gap.vectors;
  if  bools[i] = false then
      return A;
  fi;
 od;
 return true;
end;
 





