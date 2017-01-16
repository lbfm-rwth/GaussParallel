# r - rank of I_r
fourIdentityMatrices := function( r, opt... )
  local q, M;

    if not IsEmpty( opt ) then
        q := opt[1];
    else
        q := 5;
    fi;

    M := MutableCopyMat( One( GL( 2*r, q ) ) );
    M{[1..r]}{[1..r]} := One( GL( r, q ) );
    M{[r+1..2*r]}{[1..r]} := One( GL( r, q ) );
    M{[1..r]}{[r+1..2*r]} := One( GL( r, q ) );
    M{[r+1..2*r]}{[r+1..2*r]} := One( GL( r, q ) );

    #  n := MutableCopyMat( One( GL(r, q) ) );

    #randomMat := RandomInvertibleMat(r, GF(q));
    return( M );
end;
