# This file contains some help functions making the parallel parts of the code 
# in main_full_par_trafo.g more readable. (Dependency of main_full_par_trafo.g.)

# Function in step 1
ClearDownParameters := function(i, j, matrixC, TaskListClearDown, TaskListUpdateR)
	local parameters, C, D, i; #list containing galoisField, C, D, i; as defined in subprograms.g

	if (i = 1) and (j = 1) then
		C := matrixC
		
