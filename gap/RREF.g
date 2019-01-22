# Tests whether a matrix is in RREF by following the definition from
# http://mathonline.wikidot.com/reduced-row-echelon-form-of-a-matrix-rref
# According to that a matrix is in reduced row echelon form if
# 1. All of the rows that do not consist entirely of zeroes will have their first nonzero entries be 1 which we defined as leading 1s.
# 2. For any two rows that are not entirely comprised of zeroes, the leading 1 in the row below occurs farther to the right than the leading 1 in the higher rows.
# 3. Any rows consisting entirely of zeroes are placed at the bottom of the matrix.
# 4, Every column that contains a leading 1 must have zeros everywhere else in that column.

IsMatrixInRREF := function(A)
    local dimensions, i, j, k, leadingElementFound;

    dimensions := DimensionsMat(A);

    for i in [ 1 .. dimensions[1] ] do
        leadingElementFound := false;
        for j in [ 1 .. dimensions[2] ] do
            if not leadingElementFound and i < dimensions[1] then
                if not IsZero(A[i+1][j]) then
                    Info(InfoGauss, 4, "Leading ones don't go from left to right when going down");
                    return false;
                fi;
            fi;
            if not leadingElementFound and not IsZero(A[i][j]) then
                leadingElementFound := true;
                if not IsOne(A[i][j]) then
                    Info(InfoGauss, 4, "Row with a leading element that is not 1");
                    return false;
                fi;
                for k in [ 1 .. dimensions[1] ] do
                    if not k = i and not IsZero(A[k][j]) then
                        Info(InfoGauss, 4, "Column with leading one contains other elements unequal to zero");
                        return false;
                    fi;
                od;
            fi;
        od;
    od;

    return true;
end;
