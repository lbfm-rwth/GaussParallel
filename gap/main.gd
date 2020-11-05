# Declare main global functions

#! @Chapter Gaussian Elimination
#! @Section Gaussian Elimination
#!
#! @Arguments mat
#! @Returns a record
#! @Description This record contains
#! * vectors:
#!   the reduced row echelon form of the matrix mat without the zero rows
#! * heads:
#!    list that contains at position i, if nonzero, the number of the row for
#!    that the pivot element is in column i.
#! * coeffs:
#!   the transformation matrix needed to obtain the RREF from mat
#! * relations:
#!   the kernel of the matrix mat if the underlying ring is a field
#!
#! Calculates the parallel Gauss algorithm. This version assumes the underlying
#! field of the matrix by using DefaultFieldOfMatrix and uses a block size that
#! has proven to lead to the fastest execution of the algorithm. It returns the
#! commonly needed informations.
DeclareGlobalFunction( "EchelonMatTransformationBlockwise" );

#! @Arguments mat
#! @Returns a record
#! @Description This record contains
#! * vectors:
#!   the reduced row echelon form of the matrix mat without the zero rows
#! * heads:
#!    list that contains at position i, if nonzero, the number of the row for
#!    that the pivot element is in column i.
#!
#! Calculates the parallel Gauss algorithm. This version assumes the underlying
#! field of the matrix by using DefaultFieldOfMatrix and uses a block size that
#! has proven to lead to the fastest execution of the algorithm. It returns
#! a minimum of information but has the least execution time.
DeclareGlobalFunction( "EchelonMatBlockwise" );

#! @Chapter Low-Level Functions
#! @Section Low-Level Functions
#!
#!
if IsHPCGAP then
     MakeReadOnlyOrImmutableObj := MakeReadOnlyObj;
else
     MakeReadOnlyOrImmutableObj := MakeImmutable;
fi;

#! @Arguments mat, options
#! @Returns a record
#! @Description This record contains
#! * vectors:
#!   the reduced row echelon form of the matrix mat without the zero rows
#! * pivotrows and pivotcols:
#!   a list where the rows resp. columns with pivot elements are marked with 1
#!   whereas the other entries are 0
#! * rank:
#!   the rank of the matrix mat
#! * heads:
#!    list that contains at position i, if nonzero, the number of the row for
#!    that the pivot element is in column i.
#! and optionally
#! * coeffs:
#!   the transformation matrix needed to obtain the RREF from mat
#! * relations:
#!   the kernel of the matrix mat if the underlying ring is a field
#! * transformation:
#!   a transformation matrix needed to obtain the RREF from mat
#!   (filled with vectors that contain only zeros)
#!
#! Calculates the parallel Gauss algorithm. This is the version with the most
#! options to specify the functionality. You can specify how the algorithm
#! works and which return values it computes using the argument options which
#! is a record that may contain the following fields:
#! * galoisField:
#!   the galois field of the matrix, e.g. GF(7)
#! * numberBlocksHeight and numberBlocksWidth:
#!   The number of vertical and horizontal blocks in which to divide
#!   the matrix during the algorithm, note that you need to specify either none
#!   or both of those variables to make it work.
#! * numberBlocks:
#!   Use this argument if you want the same number of vertical and horizontal chops.
#! * withTrafo:
#!   A boolean specifying whether or not the transformation matrix is
#!   calculated, the default value is true.
#! * isChopped: I don't know what this does!!! TODO
#! * verify:
#!   A boolean specifying whether or not the result values should be directly
#!   checked, the default value is false.
DeclareGlobalFunction( "DoEchelonMatTransformationBlockwise" );

#! @Arguments mat
#! @Returns a record that contains information of an echelonized version of mat.
#! @Description This record contains
#! * vectors:
#!   the reduced row echelon form of the matrix mat without the zero rows
#! * heads:
#!    list that contains at position i, if nonzero, the number of the row for
#!    that the pivot element is in column i.
#! * coeffs:
#!   the transformation matrix needed to obtain the RREF from mat
#! * relations:
#!   the kernel of the matrix mat if the underlying ring is a field
#!
#! Calculates the parallel Gauss algorithm. This version assumes the underlying
#! field of the matrix by using DefaultFieldOfMatrix and uses a block size that
#! has proven to lead to the fastest execution of the algorithm. It returns the
#! commonly needed informations.
DeclareGlobalFunction( "EchelonMatTransformationBlockwise" );

#! @Arguments mat
#! @Returns a record that contains information of an echelonized version of mat.
#! @Description This record contains
#! * vectors:
#!   the reduced row echelon form of the matrix mat without the zero rows
#! * heads:
#!    list that contains at position i, if nonzero, the number of the row for
#!    that the pivot element is in column i.
#!
#! Calculates the parallel Gauss algorithm. This version assumes the underlying
#! field of the matrix by using DefaultFieldOfMatrix and uses a block size that
#! has proven to lead to the fastest execution of the algorithm. It returns
#! a minimum of information but has the least execution time.
DeclareGlobalFunction( "EchelonMatBlockwise" );
