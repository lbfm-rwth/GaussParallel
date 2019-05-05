# Declare main global functions
#! @Chapter Functions
#! @Section High-Level Functions
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
DeclareGlobalFunction("DoEchelonMatTransformationBlockwise");
DeclareGlobalFunction("EchelonMatTransformationBlockwise");
DeclareGlobalFunction("EchelonMatBlockwise");
