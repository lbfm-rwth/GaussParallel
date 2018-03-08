# Collection of small basic functions used in subfunctions of the algorithm


REX := function( galoisField,positionsBitstring,mat )
    local   i,
            up,
            down,
            upCount,
            downCount,
            numrows,
            row;
    if IsEmpty ( mat ) then return [ [],[] ]; fi;
    if IsEmpty ( positionsBitstring ) then
        ConvertToMatrixRepNC( mat,galoisField );
        return [ [],mat ];
    fi;
    numrows := Length( positionsBitstring );
    upCount := 1; 
    downCount := 1;
    up := []; 
    down := [];
    for i in [ 1 .. numrows ] do
        row := mat[ i ];
        if positionsBitstring[i] = 1 then
            up[ upCount ] := row; upCount := upCount + 1;
        else
            down[ downCount ] := row; downCount := downCount + 1;
        fi;
    od;
    ConvertToMatrixRepNC( up,galoisField );
    ConvertToMatrixRepNC( down,galoisField );

    return [ up,down ];
end;

CEX := function( galoisField,positionsBitstring,mat )
    local transposed;
    transposed := TransposedMat( mat );
    transposed := REX( galoisField,positionsBitstring,transposed );
    transposed[ 1 ] := TransposedMat( transposed[ 1 ] );
    transposed[ 2 ] := TransposedMat( transposed[ 2 ] );
    return transposed;
end;

PVC := function ( s,t )
    #We assume that positions of t correspond to zeroes in s
    local   newBitstring,
            u,
            positionU,
            positionT,
            i;
    #Case s is empty 
    if IsEmpty(s) then
        u := [];
        for i in [ 1 .. Sum( t ) ] do
            u[ i ] := 1;
        od;
        return [ t,u ];
    fi;
    if IsEmpty(t) then
        return [ s,0*[1..(Sum(s))] ];
    fi;
    #Case otherwise
    u := []; 
    newBitstring := [];
    positionU := 1;
    positionT := 1;
    for i in [ 1..Length(s) ] do
        if s[ i ] = 1 then
            newBitstring[ i ] := 1;
            u[ positionU ] := 0;
            positionU := positionU + 1;
        else
            if t[ positionT ] = 1 then
                newBitstring[ i ] := 1;
                u[ positionU ] := 1;
                positionU := positionU + 1;
            else 
                newBitstring[ i ] := 0;
            fi;
            positionT := positionT + 1;
        fi;
    od;
    return [ newBitstring,u ]; 
end;

RRF := function( galoisField,rows0,rows1,u )
    local   l,
            dim,
            index,
            index0,
            index1,
            new;
    if IsEmpty(rows0) then 
        ConvertToMatrixRepNC( rows1,galoisField );
        return rows1; 
    fi;
    if IsEmpty(rows1) then 
        ConvertToMatrixRepNC( rows0,galoisField );
        return rows0; 
    fi;
        #dim := [ DimensionsMat(rows0)[1],DimensionsMat(rows1)[1] ];
    l := Length(u);
    index0 := 1;
    index1 := 1;
    new := [];
    index := 0;
    while ( index0 + index1 -1 <= l ) do
        index := index0 + index1 -1;
        if u[index] = 0 then
            new[index] := rows0[index0];
            index0 := index0 + 1;
        else 
            new[index] := rows1[index1];
            index1 := index1 + 1;
        fi;
    od;

    ConvertToMatrixRepNC( new,galoisField );
    return new;    
end;

CRZ := function( galoisField,mat,u,nr ) 
#####################################
## !! Mind order of inputs to RRF !!
#####################################
    local   nullMat,
            numZero,
            tmp,
            sum;
    if IsEmpty(u) then return mat; fi;
    sum := Sum(u);
    if IsEmpty(mat) then
        if sum = 0 then
            return [];
        else
            return NullMat( nr,sum,galoisField ); 
        fi;
    fi;
    #numZero := Sum(u);
    ## or is it Length-Sum ? -> order of args in RRF !!
    nullMat := NullMat( sum,DimensionsMat(mat)[1],galoisField );
    return TransposedMat( RRF( galoisField,TransposedMat( mat ),nullMat,u ) );
end;

ADI := function( galoisField,mat,bitstring  )
    local   one,
            posOfNewOnes,
            i,
            copy,
            l;
    if IsEmpty(mat) then
        return IdentityMat( Length(bitstring),galoisField );
    else
        copy := MutableCopyMat( mat );
    fi;
    l := Length( bitstring );
    one := One( galoisField );
    posOfNewOnes := [];
    for i in [ 1 .. l ] do
        if bitstring[i] = 0 then
            continue;
        fi;
        Add( posOfNewOnes,i );
    od;
    l := Length( posOfNewOnes );
    for i in [ 1 .. l ] do
        copy[ i ][ posOfNewOnes[i] ] := one;
    od;
    ConvertToMatrixRepNC( copy,galoisField );
    return copy;
end;

MKR := function( bitstring,subBitstring )
    local   newBitstring,
            current,
            i,
            l;
    if IsEmpty(subBitstring) then
        return 0*[1..Sum(bitstring)]+1;
    fi;
    l := Length( bitstring  );
    newBitstring := [];
    current := 1;
    for i in [ 1 .. l ] do
        if bitstring[ i ] = 0 then
            continue;
        else
            if subBitstring[ i ] = 1 then
                newBitstring[ current ] := 0;
            else
                newBitstring[ current ] := 1;
            fi;
        current := current + 1;
        fi;
    od;
    return newBitstring;
end;

ECH := function( f,H )
    local   sct,
            Mct,
            Kct,
            tct,
            Rct,
            EMT,
            m,
            k,
            M,
            K,
            R,
            S,
            N,
            r,
            s,
            t,
            i,
            ind,
            Id,
            one,
            zero,
            dims,
            dimId;

    if H = [] then return [ [],[],[],[],[] ]; fi;
    #if Rank(H) = 0 then return [ [],[],[],[],[] ]; fi;   #--noetig?
 
    EMT := EchelonMatTransformation( H );
    m := TransposedMat(EMT.coeffs);
    if IsEmpty(m) then return [ [],[],[],[],[] ]; fi;
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
            if not IsEmpty(k) then
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

    M := -TransposedMat(M);
    ConvertToMatrixRepNC( M,f );
    K := TransposedMat(K);
    ConvertToMatrixRepNC( K,f );
    R := -TransposedMat(R);
    ConvertToMatrixRepNC( R,f );
    return [M,K,R,s,t];
 end;
