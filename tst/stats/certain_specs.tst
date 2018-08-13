gap> #There was a bug with this specification so it should be tested explicitely.
gap> isParallel := true;; dimension :=5;; rank := 1;; ring := GF(2);; numberChops := 1;;
gap> CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberChops, numberChops);;
