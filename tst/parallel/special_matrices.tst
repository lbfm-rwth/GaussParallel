gap> START_TEST("parallel/special_matrices.tst");
gap> ReadPackage("GaussPar", "tst/testdata/matrices.g");;
gap> ReadPackage("GaussPar", "tst/testfunctions.g");;
gap> for i in [1..Length(M)] do
>     result := GAUSS_TestEchelonMatTransformationBlockwiseWithGivenEchelonForm(
>     M[i], M_height[i], M_width[i], randomSource, M_q[i],
>     M_numberBlocks_height[i], M_numberBlocks_width[i], true);
>     if not result then
>         Print("Error: Special matrix number ", i);
>     fi;
> od;

# Error: No matrix
gap> DoEchelonMatTransformationBlockwise(3, rec( galoisField := GF(2), numberBlocksHeight := 2, numberBlocksWidth := 2));
Error, <mat> is not a matrix.

# Error: more blocks than rows or columns
gap> DoEchelonMatTransformationBlockwise(M1, rec( numberBlocksHeight := 3, numberBlocksWidth := 2));
Error, <numberBlocksHeight> and <numberBlocksWidth> must be less or equal than\
 the number of rows and columns respectively
gap> DoEchelonMatTransformationBlockwise(M1, rec( numberBlocksHeight := 2, numberBlocksWidth := 3));
Error, <numberBlocksHeight> and <numberBlocksWidth> must be less or equal than\
 the number of rows and columns respectively
gap> STOP_TEST("parallel/special_matrices.tst");
