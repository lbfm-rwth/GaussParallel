gap> START_TEST("stats/certain_specs.tst");
gap> #There was a bug with this specification so it should be tested explicitely.
gap> isParallel := true;; dimension :=5;; rank := 1;; ring := GF(2);; numberBlocks := 1;;
gap> GAUSS_CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberBlocks, numberBlocks);;
gap> STOP_TEST("stats/certain_specs.tst");
