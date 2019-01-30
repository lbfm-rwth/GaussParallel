gap> A := NullMat(3, 3, GF(19));;
gap> result := EchelonMatBlockwise(A);;
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> EchelonMatTransformationBlockwise(A);;
gap> B := IdentityMat(30, GF(5));;
gap> result := EchelonMatBlockwise(B);;
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> result := EchelonMatTransformationBlockwise(B);;
gap> "coeffs" in Set( RecNames( result ) );
true
gap> "relations" in Set( RecNames( result ) );
true
gap> randomSource := RandomSource(IsMersenneTwister);;
gap> C := RandomEchelonMat(50, 50, 30, randomSource, GF(13));;
gap> result := EchelonMatBlockwise(C);;
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> result := EchelonMatTransformationBlockwise(C);;
gap> "coeffs" in Set( RecNames( result ) );
true
gap> "relations" in Set( RecNames( result ) );
true
