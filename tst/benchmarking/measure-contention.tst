gap> START_TEST("stats/measure-contention.tst");
gap> n := 400;; numberBlocks := 8;; q := 5;;
gap> A := RandomMat(n, n, GF(q));;
gap> GAUSS_MeasureContention(numberBlocks, q, A, false);
gap> STOP_TEST("stats/measure-contention.tst");
