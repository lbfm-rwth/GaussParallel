# Declare global functions of this package
if IsHPCGAP then
     MakeReadOnlyOrImmutableObj := MakeReadOnlyObj;
else
     MakeReadOnlyOrImmutableObj := MakeImmutable;
fi;

DeclareGlobalFunction( "DoEchelonMatTransformationBlockwise", "Calculates echelon form and all other information" );
DeclareGlobalFunction( "EchelonMatTransformationBlockwise", "Calculates echelon form and transformation matrix of matrix" );
DeclareGlobalFunction( "EchelonMatBlockwise", "Calculates echelon form of matrix" );
