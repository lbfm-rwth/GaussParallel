# Declare global functions
#! @Chapter Gaussian Elimination
#! @Section Utility Functions
#!
#! @Arguments height, width, rank, randomSeed, ring
#! @Returns mat
#! @Description Constructs a random matrix in echelon form of the given
#! height, width and rank over the given ring. One random matrix can be
#! recreated by using the same random seed.
DeclareGlobalFunction( "RandomEchelonMat" );
#! @BeginExampleSession
#! gap> # We use Mersenne twister as a random seed here.
#! gap> randomSeed := RandomSource(IsMersenneTwister, 42);;
#! gap> M := RandomEchelonMat( 10, 10, 5, randomSeed, GF( 5 ) );;
#! gap> Display( M );
#! 1 . . . . . 1 4 4 .
#! . 1 . . . 1 3 2 4 3
#! . . 1 . . . 3 . 1 .
#! . . . 1 . 2 2 3 4 2
#! . . . . 1 2 2 4 1 4
#! . . . . . . . . . .
#! . . . . . . . . . .
#! . . . . . . . . . .
#! . . . . . . . . . .
#! . . . . . . . . . .
#! @EndExampleSession
