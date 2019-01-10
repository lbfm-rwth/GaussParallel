# Matrices in RREF
A := [[1, 0, 1, 0], [0, 1, 0, 2]];;
B := [[1, 0, 0], [0, 1, 0], [0, 0, 1]];;
C := [[1, 0], [0, 1], [0, 0], [0, 0], [0, 0]];;

# One of the leading elements is not 1
D := [[3, 0, 0, 0], [0, 1, 3, 0], [0, 0, 0, 1], [0, 0, 0, 0]];;
F := [[1, 0, 0], [0, 2, 3], [0, 0, 0]];;

# Leading zeros don't go from left to right, up to down
G := [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 0]];;
H := [[0, 1, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]];;

# Rows that contain only zeros are placed at the bottom of the matrix
J := [[0, 0, 0], [1, 4, 3], [0, 0, 0]];;
K := [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0]];;

# Column with leading 1 doesn't contain only zeros
L := [[1, 2, 0, 1, 0, 0], [0, 0, 1, 3, 1, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1]];;

matrices_rref := [ A, B, C ];;
matrices_not_rref := [ D, F, G, H, J, K, L ];;
