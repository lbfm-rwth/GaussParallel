gap> isParallel := false;;
gap> dimension := 100;;
gap> rank := 7;;
gap> ring := GF(5);;
gap> numberBlocks := 5;;
gap> average := GAUSS_CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberBlocks, numberBlocks)[3];;
gap> isParallel := true;;
gap> average := GAUSS_CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberBlocks, numberBlocks)[3];;
