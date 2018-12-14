# Matrices to be used in tests for the main algorithm.

randomSource := RandomSource(IsMersenneTwister);;

# zero matrix, small
M1 := [[0, 0], [0, 0]] * IdentityMat(2, GF(7));;

# zero matrix, large
M2 := RandomEchelonMat(100, 100, 0, randomSource, GF(13));;

# small matrix, full rank
M3 := RandomEchelonMat(5, 5, 5, randomSource, GF(5));;

# small matrix, small rank
M4 := RandomEchelonMat(5, 5, 2, randomSource, GF(11));;

# large matrix, full rank
M5 := RandomEchelonMat(200, 200, 200, randomSource, GF(17));;

# large matrix, small rank
M6 := RandomEchelonMat(200, 200, 5, randomSource, GF(3));;

# 1x1, zero rank
M7 := RandomEchelonMat(1, 1, 0, randomSource, GF(11));;

# 1x1, full rank
M8 := RandomEchelonMat(1, 1, 1, randomSource, GF(11));;

# The following matrices aren't working right now.
# That is the reason why they aren't used for testing yet.
# 1x100
M9 := RandomEchelonMat(1, 100, 1, randomSource, GF(5));

# 100x1
M10 := RandomEchelonMat(100, 1, 1, randomSource, GF(19));

# 200x50
M11 := RandomEchelonMat(200, 50, 34, randomSource, GF(7));

M_height := [2, 100, 5, 5, 200, 200, 1, 1, 1, 100, 200];;
M_width := [2, 100, 5, 5, 200, 200, 1, 1, 100, 1, 50];;
M_numberChops_height := [1, 5, 1, 1, 10, 4, 1, 1, 1, 100, 8];;
M_numberChops_width := [1, 5, 1, 1, 10, 4, 1, 1, 100, 1, 2];;
M_q := [7, 13, 5, 11, 17, 3, 11, 11, 5, 19, 7];;

M := [M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11];;
