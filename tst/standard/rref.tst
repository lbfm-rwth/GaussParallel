gap> ReadPackage("GaussPar", "tst/testdata/rref.g");;
gap> for i in [ 1 .. 3 ] do if not IsMatrixInRREF(matrices_rref[i]) then Print("Error: Matrix nr. ", i, " in RREF recognized as not in RREF\n"); fi; od;
gap> for i in [ 1 .. 7 ] do if IsMatrixInRREF(matrices_not_rref[i]) then Print("Error: Matrix nr. ", i, " not in RREF recognized as being in RREF\n"); fi; od;
