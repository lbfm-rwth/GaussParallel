gap> START_TEST("standard/echelonbasismutable.tst");

# Some random full rank matrices. Note that a random mat is usually full rank.
gap> for i in [1..20] do
> A := CVEC_RandomMat(50, 50, 5, 2);;
> basis := GAUSS_EchelonBasisMutableT(A);;
> sortPerm := SortingPerm(basis!.pivots);;
> m := Permuted(CMat(basis!.vectors), sortPerm);;
> if IsOne(m) then Print(i, ", "); fi;
> od; Print("\n");
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 

# Corner cases
gap> zeros := [];;
gap> Add(zeros, CVEC_ZeroMat(0, 0, 5, 2));
gap> Add(zeros, CVEC_ZeroMat(1, 1, 5, 2));
gap> Add(zeros, CVEC_ZeroMat(0, 5, 5, 2));
gap> Add(zeros, CVEC_ZeroMat(5, 0, 5, 2));
gap> Add(zeros, CVEC_ZeroMat(5, 5, 5, 2));
gap> Perform(zeros, GAUSS_EchelonBasisMutableTX);
gap> ids := [];;
gap> Add(ids, CVEC_IdentityMat(0, 5, 2));
gap> Add(ids, CVEC_IdentityMat(1, 5, 2));
gap> Add(ids, CVEC_IdentityMat(5, 5, 2));
gap> Perform(ids, GAUSS_EchelonBasisMutableTX);

# Low rank matrices
gap> for i in [1..20] do
> m := CVEC_ZeroMat(100, 100, 5, 3);;
> CopySubMatrix(CVEC_IdentityMat(30, 5, 3), m,
>               [1..30], [1..30], [1..30], [1..30]);
> if Length(Vectors(GAUSS_EchelonBasisMutableT(m))) = RankMat(m) then
>    Print(i, ", ");
> fi;
> od; Print("\n");
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 

# Permutation matrices
gap> for i in [1..20] do
> perm := PseudoRandom(SymmetricGroup(100));;
> permMat := Permuted(CVEC_IdentityMat(100, 5, 1), perm);;
> res := GAUSS_EchelonBasisMutableT(permMat);
> if perm = SortingPerm(res!.pivots) ^ -1 then Print(i, ", "); fi;
> od; Print("\n");
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 

# STOP
gap> STOP_TEST("standard/echelonbasismutable.tst");
