# Declare global functions
#! @Chapter Functions
#! @Section Low-Level Functions
#!
#! @Arguments mat
#! @Returns bool
#! @Description Checks whether or not a matrix mat is in RREF
#! where RREF means 'Reduced Row Echelon Form'. A matrix is in RREF if it
#! meets the following criteria:
#! * The first nonzero entry of a nonzero row is 1 (called 'leading 1').
#! * When looking at two nonzero rows, than the leading one of the lower
#!   rows is further right than the other leading one.
#! * The rows containing only zeros are below the nonzero rows.
#! * Every column that contains a leading one must have zeros everywhere else
#!   in that column.
DeclareGlobalFunction( "IsMatrixInRREF" );
#! @BeginExampleSession
#! gap> M := RandomMat( 3, 3 );;
#! gap> Display( M );
#! [ [   1,   0,  -1 ],
#!   [  -1,  -1,   1 ],
#!   [  -1,   1,  -2 ] ]
#! gap> IsMatrixInRREF( M );
#! false
#! gap> N := IdentityMat( 3 );;
#! gap> Display( N );
#! [ [  1,  0,  0 ]
#!   [  0,  1,  0 ],
#!   [  0,  0,  1 ] ]
#! gap> IsMatrixInRREF( N );
#! true
#! gap> L := [ [ 1, 0, 3 ], [ 0, 1, 5 ] ];;
#! gap> Display( L );
#! [ [  1,  0,  3 ],
#!   [  0,  1,  5 ] ]
#! gap> IsMatrixInRREF( L );
#! true
#! @EndExampleSession
