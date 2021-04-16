# Declare global functions
#! @Chapter Utility Functions
#! @Section Utility Functions
#!
#! @Arguments mat
#! @Returns bool
#! @Description Checks whether the matrix **mat** is in RREF
DeclareGlobalFunction( "IsMatrixInRREF" );
#! @BeginExampleSession
#! gap> M := RandomMat( 3, 3 );;
#! gap> Display( M );
#! [ [   1,   0,  -1 ],
#!   [  -1,  -1,   1 ],
#!   [  -1,   1,  -2 ] ]
#! gap> IsMatrixInRREF( M );
#! false
#! gap> L := [ [ 1, 0, 3 ], [ 0, 1, 5 ] ];;
#! gap> Display( L );
#! [ [  1,  0,  3 ],
#!   [  0,  1,  5 ] ]
#! gap> IsMatrixInRREF( L );
#! true
#! @EndExampleSession
