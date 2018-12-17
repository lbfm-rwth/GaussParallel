gap> n := 400;; numberBlocks := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> MeasureContention(numberBlocks, q, A, false);
