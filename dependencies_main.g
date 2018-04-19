# This file contains some help functions making the parallel parts of the code 
# in main_full_par_trafo.g more readable. (Dependency of main_full_par_trafo.g.)

# Function in step 1
ClearDownParameters := function(i, j, matrixC, TaskListClearDown, TaskListUpdateR, galoisField, dummyTask)
	local list, C, D, number; # parameters for ClearDown as in subprograms.g

	if (i = 1) and (j = 1) then
		list := [ dummyTask ];
		C := matrixC[1][1];
		D := rec( vectors:=[], bitstring:=[] );
	elif (i = 1) and (j > 1) then
		list := [ TaskListUpdateR[1][j-1][j] ];
		C := TaskResult( TaskListUpdateR[1][j-1][j] ).C;
		D := rec( vectors:=[], bitstring:=[] );
	elif (i > 1) and (j = 1) then
		list := [ TaskListClearDown[i-1][1] ];
		C := matrixC[i][1];
		D := TaskResult( TaskListClearDown[i-1][1] ).D;
	else
		list := [ TaskListClearDown[i-1][j], TaskListUpdateR[i][j-1][j] ];
		C := TaskResult( TaskListUpdateR[i][j-1][j] ).C;
		D := TaskResult( TaskListClearDown[i-1][j] ).D;
	fi;
		
	number := i;
	return [ list, galoisField, C, D, number ];
end;

ExtendParameters := function(i, j, TaskListClearDown, TaskListE)
	local list, A, E, flag; # parameters for Extend as in subprograms.g
	
	if (i = 1) and (j = 1) then
		list := [ TaskListClearDown[i][j] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		E := rec( rho := [], delta := [], nr := 0 );
	elif (i = 1) and (j > 1) then
		list := [ TaskListClearDown[i][j], TaskListE[i][j-1] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		E := TaskResult( TaskListE[i][j-1] );
	elif (i > 1) and (j = 1) then
		list := [ TaskListClearDown[i][j] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		E := rec( rho := [], delta := [], nr := 0 );
	else
		list := [ TaskListClearDown[i][j], TaskListE[i][j-1] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		E := TaskResult( TaskListE[i][j-1]);
	fi;

	flag := j;
	return [ list, A, E, flag ];
end;

UpdateRowParameters := function(i, j, k, matrixC, TaskListClearDown, TaskListUpdateR, galoisField)
	local list, A, C, B, number; # parameters for UpdateRow as in subprograms.g

	if (i = 1) and (j = 1) then
		list := [ TaskListClearDown[i][j] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		C := matrixC[i][k];
		B := [];
	elif (i = 1) and (j > 1) then
		list := [ TaskListClearDown[i][j], TaskListUpdateR[i][j-1][k] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		C := TaskResult( TaskListUpdateR[i][j-1][k] ).C;
		B := [];
	elif (i > 1) and (j = 1) then
		list := [ TaskListClearDown[i][j], TaskListUpdateR[i-1][j][k] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		C := matrixC[i][k];
		B := TaskResult( TaskListUpdateR[i-1][j][k] ).B;
	else
		list := [ TaskListClearDown[i][j], TaskListUpdateR[i-1][j][k], TaskListUpdateR[i][j-1][k] ];
		A := TaskResult( TaskListClearDown[i][j] ).A;
		C := TaskResult( TaskListUpdateR[i][j-1][k] ).C;
		B := TaskResult( TaskListUpdateR[i-1][j][k] ).B;
	fi;

	number := i;
	return [ list, galoisField, A, C, B, number ];
end;
