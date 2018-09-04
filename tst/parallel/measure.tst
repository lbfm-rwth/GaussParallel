gap> n := 400;; numberChops := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> MeasureContention(numberChops, q, A, false);
