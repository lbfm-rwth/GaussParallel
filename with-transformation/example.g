# Arbitrary 20 x 30 Matrix over GF(5), not having full rank, containing nullrows and collumns

ExampleGaussParallel := function(  )
 local gap,A,l;

 A := RandomMat(10,30,GF(5));
 A := Concatenation( NullMat(5,30,GF(5)),A );
 A := Concatenation( A,RandomMat( 5,30,GF(5) ) );
 A [12] := A[12]*Zero(GF(5));
 A [13] := A[12]*Zero(GF(5));
 A := MutableCopyMat(TransposedMat(A));
 A[17] := A[17]*Zero(GF(5));
 A[18] := A[17]*Zero(GF(5));
 A[19] := A[17]*Zero(GF(5));
 A := TransposedMat(A);

 Print("Echolonizing the following matrix over GF(5):");
 Print("\n");
 Display(A);
 
 l := GaussParallel( A );
 gap := EchelonMatTransformation( A );

 Print("Echolonized matrix:");
 Print("\n");
 Display(l.vectors);

 Print("Transformation of pivotrows:");
 Print("\n");
 Display(l.coeffs);
 
 Print("Transformation of non-pivotrows to zero:");
 Print("\n");
 Display(l.relations);

 Print("Results agree with gap: ");
 Print(l.vectors = gap.vectors and l.coeffs = gap.coeffs and l.relations = gap.relations);

 Print("\n");
 return 0;
end;
