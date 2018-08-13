gap> isParallel := false;;
gap> dimension := 100;;
gap> rank := 7;;
gap> ring := GF(5);;
gap> numberChops := 5;;
gap> average := CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberChops, numberChops)[3];;
gap> isParallel := true;;
gap> average := CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberChops, numberChops)[3];;
