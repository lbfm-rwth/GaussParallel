if IsHPCGAP then
     MakeReadOnlyOrImmutableObj := MakeReadOnlyObj;
else
     MakeReadOnlyOrImmutableObj := MakeImmutable;
fi;

# Declare main global functions
#! @Chapter Gaussian Elimination
#! @Section Gaussian Elimination


#! @Arguments mat, [options]
#! @Returns a record that contains information on the echelon form of **mat** and the corresponding transformation matrix.
#! @Description This is the main function of the GaussPar package. It computes the reduced row echelon form (RREF) of the matrix **mat** and the corresponding transformation matrix. In a pre-processing step, **mat** is split up into a block matrix whose blocks can be processed in parallel.
#!
#!  The input parameters have the following meaning:
#! * **mat** is a matrix defined over a finite field
#! * **options** is a record that can be used to provide some additional parameters.  Note that a specification of either **numberBlocks** or both of
#!   **numberBlocksHeight** and **numberBlocksWidth** is mandatory. The following parameters are currently supported:
#!   * **numberBlocksHeight** and **numberBlocksWidth**:
#!   The number of vertical and horizontal blocks in which to divide
#!   the matrix during the algorithm.
#!   * **numberBlocks**:
#!   Use this argument if you want the same number of vertical and horizontal blocks in the block matrix decomposition of **mat**.
#!   * **verify**: If set to **true**, the computation is verified at the end. That is, we check wheter **coeffs** * **mat** is in RREF. This option is only available for the function **EchelonMatTransformationBlockwise**.
#!   * **isChopped**: It is possible to input **mat** directly as a matrix of block matrices. In this case the parameter **isChopped** must be set to **true** and the splitting step is skipped.
#!
#!  The output record contains the following items:
#! * **vectors**:
#!   a matrix that forms the RREF of **mat** without zero rows
#! * **heads**:
#!    a list of integers, such that **heads[i]** gives the number of the row for
#!    which the pivot element is in column i. If no such row exists, **heads[i]** is **0**.
#! * **coeffs**:
#!   the corresponding transformation matrix. It holds **coeffs** * **mat** = **vectors**.
#! * **relations**: a matrix whose rows form a basis for the row null space of **mat**. If **relations** is not the  empty list, it holds that **relations** * **mat** = **0**. Otherwise **mat** has full row rank.
#!  
#!
#! @BeginExampleSession
#! gap> A := RandomMat(8, 5, GF(5)) * RandomMat(5, 8, GF(5));;
#! gap> Display(A);
#! 1 4 3 2 4 4 3 4
#! 4 1 2 4 2 . . 4
#! 2 3 1 4 3 3 1 3
#! 3 . 4 3 3 2 4 .
#! 4 1 3 2 3 3 . 2
#! 2 1 3 3 1 1 2 3
#! . 3 3 . 1 1 3 .
#! 4 1 4 1 4 3 1 1
#! gap> res := EchelonMatTransformationBlockwise(A, rec(numberBlocks := 2));;
#! gap> Display(res.vectors);
#! 1 . . . . . . 3
#! . 1 . . . 3 2 3
#! . . 1 . . . . 3
#! . . . 1 . 2 1 1
#! . . . . 1 2 2 2
#! gap> res.coeffs * A=res.vectors;
#! true
#! @EndExampleSession
#!  The transformation matrix can be easily obtained from the output record as follows:
#! @BeginExampleSession
#! gap> trafo := Concatenation(res.coeffs, res.relations);;
#! gap> Display(trafo * A);
#! 1 . . . . . . 3
#! . 1 . . . 3 2 3
#! . . 1 . . . . 3
#! . . . 1 . 2 1 1
#! . . . . 1 2 2 2
#! . . . . . . . .
#! . . . . . . . .
#! . . . . . . . .
#! @EndExampleSession
#!It is also possible to directly input a block matrix.
#! @BeginExampleSession
#! gap> A := RandomMat(8, 5, GF(5)) * RandomMat(5, 8, GF(5));;
#! gap> res1 := EchelonMatTransformationBlockwise(A, rec(numberBlocks := 2));;
#! gap> A1 := A{[1..4]}{[1..4]};;
#! gap> A2 := A{[1..4]}{[5..8]};;
#! gap> A3 := A{[5..8]}{[1..2]};;
#! gap> A4 := A{[5..8]}{[5..8]};;
#! gap> A_blockwise := [[A1,A2],[A3,A4]];;
#! gap> res3 := EchelonMatTransformationBlockwise(A_blockwise, 
#!   rec(numberBlocks := 2, isChopped := true));;
#! gap> res3 = res1;
#! true
#! @EndExampleSession
DeclareGlobalFunction( "EchelonMatTransformationBlockwise" );

#! @Arguments mat, [options]
#! @Returns a record that contains information on the echelon form of **mat**.
#! @Description This is a version of the main function that computes the reduced row echelon form (RREF) of the matrix **mat** but doesn't compute the corresponding transformation matrix. In a pre-processing step, **mat** is split up into a block matrix whose blocks can be processed in parallel.
#! 
#!  The input parameters have the following meaning:
#! * **mat** is a matrix defined over a finite field
#! * **options** is a record that is used to provide additional parameters. For more information see **EchelonMatTransformationBlockwise** (See  <Ref Func="EchelonMatTransformationBlockwise"/>).
#!  
#!  The output record contains the items **vectors** and **heads**. For their meaning, see <Ref Func="EchelonMatTransformationBlockwise"/>.
DeclareGlobalFunction( "EchelonMatBlockwise" );

DeclareGlobalFunction( "DoEchelonMatTransformationBlockwise" );
