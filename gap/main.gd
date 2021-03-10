# Declare main global functions

#! @Chapter Gaussian Elimination
#! @Section Gaussian Elimination


#! @Arguments mat, [options]
#! @Returns a record that contains information on the echelon form of **mat** and the corresponding transformation matrix.
#! @Description This is the main function of the GaussPar package. It computes the reduced row echelon form (RREF) of the matrix **mat** and the corresponding transformation matrix. In a pre-processing step, **mat** is split up into smaller block matrices which can be processed in parallel.
#!  
#!  The output record contains the following items:
#! * **vectors**:
#!   a matrix that forms the RREF of **mat** without zero rows
#! * **heads**:
#!    a list of integers, such that **heads[i]** gives the number of the row for
#!    which the pivot element is in column i. If no such row exists, **heads[i]** is **0**.
#! * **coeffs**:
#!   the corresponding transformation matrix. It holds **coeffs** * **mat** = **vectors**.
#! * **relations**: the kernel of the matrix **mat**. If **relations** is not the  empty list, it holds **relations** * **mat** = **0**. Otherwise **mat** has full row rank.
#! 
#!  The input parameters have the following meaning:
#! * **mat** is a matrix defined over a finite field
#! * **options** is a record that can be used to provide some additional parameters. The following   are currently supported:
#! * **numberBlocksHeight** and **numberBlocksWidth**:
#!   The number of vertical and horizontal blocks in which to divide
#!   the matrix during the algorithm, note that you need to specify either none
#!   or both of those variables to make it work.
#! * **numberBlocks**:
#!   Use this argument if you want the same number of vertical and horizontal chops.
#!   * **verify**: If set to **true**, the computation will be verified at the end.
#!   * **isChopped**: It is possible to input **mat** directly as a matrix of block matrices. In this case the parameter **isChopped** must be set to **true** and the splitting step is skipped. Note that a specification of **numberBlocks** is mandatory if **isChopped** is set to **true**. 
DeclareGlobalFunction( "EchelonMatTransformationBlockwise" );

#! @Arguments mat, [options]
#! @Returns a record that contains information on the echelon form of **mat**.
#! @Description This is a version of the main function that computes the reduced row echelon form (RREF) of the matrix **mat** but doesn't compute the corresponding transformation matrix. In a pre-processing step, **mat** is split up into smaller block matrices which can be processed in parallel.
#!  
#!  The output record contains the following items:
#! * **vectors**:
#!   a matrix that forms the RREF of **mat** without zero rows
#! * **heads**:
#!    a list of integers, such that **heads[i]** gives the number of the row for
#!    which the pivot element is in column i. If no such row exists, **heads[i]** is **0**.
#! 
#!  The input parameters have the following meaning:
#! * **mat** is a matrix defined over a finite field
#! * **options** is a record that can be used to provide some additional parameters. For more information see **EchelonMatTransformationBlockwise** (TODO:Use cross-ref here).
DeclareGlobalFunction( "EchelonMatBlockwise" );

#! @Chapter Low-Level Functions
#! @Section Low-Level Functions
#!
#!
#! @Arguments mat, options
#! @Returns a record that contains information on the echelon form of **mat**.
#! @Description This record contains
#! * **vectors**:
#!   the reduced row echelon form of the matrix **mat** without the zero rows
#! * **pivotrows** and **pivotcols**:
#!   a list where the rows resp. columns with pivot elements are marked with 1
#!   whereas the other entries are 0
#! * **rank**:
#!   the rank of the matrix **mat**
#! * **heads**:
#!    list that contains at position i, if nonzero, the number of the row for
#!    that the pivot element is in column i.
#! and optionally
#! * **coeffs**:
#!   the transformation matrix needed to obtain the RREF from mat
#! * **relations**:
#!   the kernel of the matrix mat if the underlying ring is a field
#! * **transformation**:
#!   a transformation matrix needed to obtain the RREF from mat
#!   (filled with vectors that contain only zeros)
#!
#! This function executes the parallel Gaussian algorithm. Both **EchelonMatTransformationBlockwise** and **EchelonMatBlockwise** are wrapper functions of it.
#! This is the version with the most
#! options to specify functionality. You can specify how the algorithm
#! works and which return values it computes using the argument options which
#! is a record that may contain the following items:
#! * **galoisField**:
#!   the galois field of the matrix, e.g. GF(7)
#! * **numberBlocksHeight** and **numberBlocksWidth**:
#!   The number of vertical and horizontal blocks in which to divide
#!   the matrix during the algorithm, note that you need to specify either none
#!   or both of those variables to make it work.
#! * **numberBlocks**:
#!   Use this argument if you want the same number of vertical and horizontal chops.
#! * **withTrafo**:
#!   A boolean specifying whether or not the transformation matrix is
#!   computed, the default value is true.
#!   * **isChopped**: It is possible to input **mat** directly as a matrix of block matrices. In this case the parameter **isChopped** must be set to **true** and the splitting step is skipped. Note that a specification of **numberBlocks** is mandatory if **isChopped** is set to **true**.
#! * **verify**:
#!   A boolean specifying whether or not the result values should be directly
#!   checked, the default value is false.
DeclareGlobalFunction( "DoEchelonMatTransformationBlockwise" );
