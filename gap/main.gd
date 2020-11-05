# Declare global functions of this package
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
DeclareGlobalFunction("DoEchelonMatTransformationBlockwise");
DeclareGlobalFunction("EchelonMatTransformationBlockwise");
DeclareGlobalFunction("EchelonMatBlockwise");
