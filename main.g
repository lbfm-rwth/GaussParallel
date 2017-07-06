### TODO ###
## rename variables?, e.g.:
## t -> ListOfPivotColumns
### TODO ###
GaussParallel := function( A )
    local i, ct, k, f, n, smallDivisors, l, rank, nrows, ncols;
    f := DefaultFieldOfMatrix( A );
    nrows := DimensionsMat(A)[1];
    ncols := DimensionsMat(A)[2];
    n := Gcd( nrows, ncols );

    # Error( "Break Point - before Step1" );
    l := Step1( A, n );
    l := Step2( f, n, nrows, ncols, l );
    return rec(vectors := -l[1], columnPermutation := l[2] );
end;
