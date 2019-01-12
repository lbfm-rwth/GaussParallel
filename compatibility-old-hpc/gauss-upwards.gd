#############################################################################
##
##  GaussDense.gd                Gauss package                Simon Goertzen
##
##  Copyright 2007-2008 Lehrstuhl B für Mathematik, RWTH Aachen
##
##  Declaration stuff for Gauss algorithms on dense (IsMatrix) matrices.
##
#############################################################################

##
DeclareOperation( "EchelonMatTransformationDestructive", #RREF over a ring, returns the same record as SemiEchelonMatTransformation but with ordered vectors
        [ IsMatrix ] );

DeclareAttribute( "EchelonMatTransformation",
        IsMatrix );
