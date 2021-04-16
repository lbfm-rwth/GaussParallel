gap> START_TEST("standard/main.tst");
gap> A := NullMat(3, 3, GF(19));;
gap> result := EchelonMatBlockwise(A, rec(numberBlocks := 1));;
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> EchelonMatTransformationBlockwise(A, rec(numberBlocks := 1));;
gap> B := IdentityMat(30, GF(5));;
gap> result := EchelonMatBlockwise(B, rec(numberBlocks := 2));;
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> result := EchelonMatTransformationBlockwise(B, rec(numberBlocks := 2));;
gap> "coeffs" in Set( RecNames( result ) );
true
gap> "relations" in Set( RecNames( result ) );
true
gap> randomSource := RandomSource(IsMersenneTwister);;
gap> C := RandomEchelonMat(50, 50, 30, randomSource, GF(13));;
gap> result := EchelonMatBlockwise(C, rec(numberBlocks := 2));;
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> result := EchelonMatTransformationBlockwise(C, rec(numberBlocks := 2));;
gap> "coeffs" in Set( RecNames( result ) );
true
gap> "relations" in Set( RecNames( result ) );
true
gap> STOP_TEST("standard/main.tst");
