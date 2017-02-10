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
    smallDivisors := Filtered( [5,4,3,2], x -> n mod x = 0 );
    if not IsEmpty( smallDivisors ) then
        n := smallDivisors[1];
    else
        Info( InfoWarning, 1, "Using gcd to chop A. May be too big!" );
    fi;

    # Error( "Break Point - before Step1" );
    l := Step1( A, n );
    l := Step2( f, n, nrows, ncols, l );
    return rec(vectors := -l[1], columnPermutation := l[2] );
end;
