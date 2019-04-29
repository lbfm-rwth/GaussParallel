# Declare main global functions
#! @Chapter High-level functions
#! @Section Test
#!
#!
if IsHPCGAP then
     MakeReadOnlyOrImmutableObj := MakeReadOnlyObj;
else
     MakeReadOnlyOrImmutableObj := MakeImmutable;
fi;
#!
#! @Description Some functions
#! @Returns Stuff
DeclareGlobalFunction( "DoEchelonMatTransformationBlockwise", "Calculates echelon form and all other information" );
DeclareGlobalFunction( "EchelonMatTransformationBlockwise", "Calculates echelon form and transformation matrix of matrix" );
DeclareGlobalFunction( "EchelonMatBlockwise", "Calculates echelon form of matrix" );
