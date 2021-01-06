# These functions are currently unused under HPCGAP because cvec causes a
# segfault. If they work then one day they should be added to the cvec package.
LoadPackage("cvec");

GAUSS_EchelonBasisMutableTX := function(m)
    local downBasis, coeffs, vectors, pivots, findPivotsPerm, nrRows, dec,
    upBasis, reversePerm, s, transformation, i;
    # Clear downwards
    downBasis := SemiEchelonBasisMutableT(m);
    if Length(m) = 0 then
        return downBasis;
    fi;
    # Unpack downBasis
    pivots := Pivots(downBasis);
    vectors := Vectors(downBasis);
    nrRows := NrRows(vectors);
    # Permuted(vectors, findPivotsPerm ^ -1) is in semi echelon form
    findPivotsPerm := SortingPerm(pivots) ^ -1;
    reversePerm := PermList(Reversed([1..nrRows]));
    # Clear upwards: Start with the row with the right-most pivot and then work
    # your way upwards.
    # We use upBasis and dec for bookkeeping.
    upBasis := EmptySemiEchelonBasis(vectors);
    dec := ZeroVector(nrRows,vectors);
    for i in [0..nrRows - 1] do
        # CleanRow updates downBasis.vectors. We update downBasis.coeffs
        # manually. downBasis.relations does not need to be updated.
        # vectors[i ^ findPivotsPerm] is the i-th pivot row.
        CleanRow(upBasis, vectors[(nrRows - i) ^ findPivotsPerm], true, dec);
        # At all times we maintain: coeffs * mat = vectors.
        # vectors is cleaned via row operations which we repeat on coeffs.
        # Let v := vectors[(nrRows - i) ^ findPivotsPerm] and upVecs :=
        # Vectors(upBasis). Then the information stored in dec is:
        # dec[i + 1] * upVecs[i + 1] =
        #    dec{[1..i]} * upVecs{[1..i]} - v.
        # Note that the downBasis.vectors are inserted into upBasis.vectors in
        # descending order wrt to pivot positions, i.e. the row with the
        # right-most pivot comes first. We have (at the end of this loop):
        # upBasis.vectors =
        #    Permuted(vectors, findPivotsPerm ^ -1 * reversePerm).
        s := -dec[i+1] ^ -1;
        dec[i+1] := -One(s);
        transformation := s * Permuted(dec, reversePerm * findPivotsPerm);
        downBasis!.coeffs[(nrRows - i) ^ findPivotsPerm] :=
            transformation * downBasis!.coeffs;
    od;
    return downBasis;
end;

GAUSS_EchelonBasisMutableT := function(m)
    local copyM;
    copyM := MutableCopyMat(m);
    return GAUSS_EchelonBasisMutableTX(copyM);
end;
