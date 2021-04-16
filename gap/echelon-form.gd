# Declare global functions
#! @Chapter Utility Functions
#! @Section Utility Functions
#!
#! @Arguments height, width, rank, randomSeed, field
#! @Returns mat
#! @Description Constructs a random matrix in RREF of the given
#! height (number of rows), width (number of columns) and rank over the given field. Using the same random seed will recreate the same matrix.
DeclareGlobalFunction( "RandomEchelonMat" );
#! @BeginExampleSession
#! gap> # We use Mersenne twister as a random seed here.
#! gap> randomSeed := RandomSource(IsMersenneTwister, 42);;
#! gap> M := RandomEchelonMat(10, 10, 5, randomSeed, GF(5));;
#! gap> Display(M);
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
