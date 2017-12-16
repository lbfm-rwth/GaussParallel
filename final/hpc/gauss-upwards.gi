##
InstallMethod( EchelonMatTransformationDestructive,
        "generic method for matrices",
        [ IsMatrix and IsMutable ],
        function( mat )
    local zero,      # zero of the ring of <mat>
          nrows,     # number of rows in <mat>
          ncols,     # number of columns in <mat>
          vectors,   # list of basis vectors
          field,
          heads,     # list of pivot positions in 'vectors'
          i,         # loop over rows
          j,         # loop over columns
          T,         # transformation matrix
          coeffs,    # list of coefficient vectors for 'vectors'
          relations, # basis vectors of the null space of 'mat'
          row,
          head,
          x,
          row2,
          rank,
          list,
	  a;
    
    nrows := Length( mat );
    ncols := Length( mat[1] );
    
    zero  := Zero( mat[1][1] );
    
    heads   := ListWithIdenticalEntries( ncols, 0 );
    vectors := [];
    
    if Characteristic( zero ) mod 2 = 0 then
        if IsGF2VectorRep( mat[1] ) then
            field := GF( 2 );
        else
            field := GF( 2^Q_VEC8BIT( mat[1] ) );
        fi;
    else
        field := Field( zero );
    fi;
    
    T         := IdentityMat( nrows, field );
    coeffs    := [];
    relations := [];
    
    for i in [ 1 .. nrows ] do
        
        row := mat[i];
        row2 := T[i];
        
        # Reduce the row with the known basis vectors.
        for j in [ 1 .. ncols ] do
            head := heads[j];
            if head <> 0 then
                x := - row[j];
                if x <> zero then
                    AddRowVector( row2, coeffs[ head ],  x );
                    AddRowVector( row,  vectors[ head ], x );
                fi;
            fi;
        od;
        
        
        j:= PositionNot( row, zero );
        if j <= ncols then
            
            # We found a new basis vector.
            x := Inverse( row[j] );
            if x = fail then
                TryNextMethod();
            fi;
            Add( coeffs,  row2 * x );
            Add( vectors, row  * x );
            heads[j]:= Length( vectors );
            
        else
            Add( relations, row2 );
        fi;
        
    od;
    
    # gauss upwards:
    
    list := Filtered( heads, x->x<>0 );
    rank := Length( list );
    
    for j in [ncols,ncols-1..1] do
        head := heads[j];
        if head <> 0 then
            a := Difference( [1..head-1], heads{[j+1..ncols]} );
            for i in a do
                row := vectors[i];
                row2 := coeffs[i];
                x := - row[j];
                if x <> zero then
                    AddRowVector( row2, coeffs[head], x );
                    AddRowVector( row, vectors[head], x );
                fi;
            od;
        fi;
    od;
    
    #order rows:
    
    vectors := vectors{list};
    
    coeffs := coeffs{list};
    
    list := Filtered( [1..ncols], j -> heads[j] <> 0 );
    heads{list} := [1..rank];  #just for compatibilty, vectors are ordered already
    
    return rec( heads := heads,
                vectors := vectors,
                coeffs := coeffs,
                relations := relations );
    
end );
##
InstallMethod( EchelonMatTransformation,
        "generic method for matrices",
        [ IsMatrix ],
        function( mat )
    local copymat, v, vc, f;
    copymat := [];
    f := DefaultFieldOfMatrix(mat);
    for v in mat do
        vc := ShallowCopy(v);
        ConvertToVectorRepNC(vc,f);
        Add(copymat, vc);
    od;
    return EchelonMatTransformationDestructive( copymat );
end);
