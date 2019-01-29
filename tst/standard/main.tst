gap> START_TEST("standard/main.tst");
gap> A := NullMat(3, 3, GF(19));;
gap> result := EchelonMatBlockwise(A);
rec( heads := [ 0, 0, 0 ], vectors := [  ] )
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> EchelonMatTransformationBlockwise(A);
rec( coeffs := [  ], heads := [ 0, 0, 0 ], 
  relations := [ [ Z(19)^0, 0*Z(19), 0*Z(19) ], [ 0*Z(19), Z(19)^0, 0*Z(19) ],
      [ 0*Z(19), 0*Z(19), Z(19)^0 ] ], vectors := [  ] )
gap> B := IdentityMat(30, GF(5));;
gap> result := EchelonMatBlockwise(B);
rec( heads := [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  vectors := < immutable compressed matrix 30x30 over GF(5) > )
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
gap> result := EchelonMatBlockwise(C);
rec( heads := [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0 ], vectors := < immutable compressed matrix 30x
    50 over GF(13) > )
gap> "coeffs" in Set( RecNames( result ) );
false
gap> "relations" in Set( RecNames( result ) );
false
gap> result := EchelonMatTransformationBlockwise(C);;
gap> "coeffs" in Set( RecNames( result ) );
true
gap> "relations" in Set( RecNames( result ) );
true
gap> STOP_TEST("standard/main.tst");
