# Matrices to be used in tests for the main algorithm.

randomSource := RandomSource(IsMersenneTwister);;

# zero matrix, small
M_zero_small := [[0, 0], [0, 0]] * IdentityMat(2, GF(7));;

# zero matrix, large
M_zero_large := RandomEchelonMat(100, 100, 0, randomSource, GF(13));;

# small matrix, full rank
M_small_full := RandomEchelonMat(5, 5, 5, randomSource, GF(5));;

# small matrix, small rank
M_small_small := RandomEchelonMat(5, 5, 2, randomSource, GF(11));;

# large matrix, full rank
M_large_full := RandomEchelonMat(200, 200, 200, randomSource, GF(17));;

# large matrix, small rank
M_large_small := RandomEchelonMat(200, 200, 5, randomSource, GF(3));;

# no matrix
M_no := 3;;

M := [M_zero_small, M_zero_large, M_small_full, M_small_small, M_large_full, M_large_small, M_no];;
M_width := [2, 100, 5, 5, 200, 200, 5];;
M_height := [2, 100, 5, 5, 200, 200, 5];;
M_numberChops := [1, 5, 1, 1, 10, 4, 1];;
M_q := [7, 13, 5, 11, 17, 3, 23];;
