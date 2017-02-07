### TODO ###
## rename variables, e.g.:
## t -> ListOfPivotColumns
### TODO ###
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

## TODO make a new test function
## This only works for sequential version (before 07 Feb 17)
#TestGaussParallel := function( nr,nc,iter )
# local i,test,A,bools,gap;
# 
# bools:=[]; 
# for i in [1..iter] do
#  A := RandomMat(3*nr,nc,GF(5));
#  A[1] := Zero(GF(5))*A[1];
#  A := RandomMat(nc,2*nr,GF(5)) * A;
#  
#  test := GaussParallel( A );
#  gap := EchelonMat(A);
#
#  bools[i] := test.vectors = gap.vectors;
#  if  bools[i] = false then
#      return A;
#  fi;
# od;
# return true;
#end;
