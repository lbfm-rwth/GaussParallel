# Declare global functions
#! @Chapter Functions
#! @Section Low-Level Functions
#!
#! @Arguments height, width, rank, randomSeed, ring
#! @Returns mat
#! @Description Constructs a random matrix in echelon form of the given
#! height, width and rank over the given ring. One random matrix can be
#! recreated by using the same random seed.
DeclareGlobalFunction( "RandomEchelonMat" );
