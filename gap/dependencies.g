# This file contains the help functions that manage the parallelization of the
# algorithm. In particular, the functions in this file manage the dependencies
# of the algorithm's subprograms between each other.

# This function is temporary until it is not necessary anymore to find
# a number of chops that divides the width or height.
# When different chop sizes are possible, we will presumably find the perfect
# chop size and divide in chops of this size.
GAUSS_calculateChops := function( a )
    # a is either height or width of the matrix
    local i;
    i := 13;

    while (i > 1) do
        if (a mod i = 0) then
            return a/i;
        fi;
        i := i - 1;
    od;

    return 1;
end;

# Function in step 1
GAUSS_ClearDownParameters := function(i, j, matrixC, TaskListClearDown, TaskListUpdateR, galoisField)
	local list, C, D, number; # parameters for GAUSS_ClearDown as in subprograms.g
    
    # The output of this function is used by ScheduleTask which has a bug
    # so that sometimes it doesn't work with an empty list. That's why
    # we currently need this workaround.
    list := [RunTask( function() return []; end )];

    if (i = 1) then
		D := rec( vectors:=[], bitstring:=[] );
    else
		D := TaskResult( TaskListClearDown[i-1][j] ).D;
		Add(list, TaskListClearDown[i-1][j]);
    fi;

    if (j = 1) then
		C := matrixC[i][1];
    else
		C := TaskResult( TaskListUpdateR[i][j-1][j] ).C;
		Add(list, TaskListUpdateR[i][j-1][j]);
    fi;
		
	number := i;
	return rec( dependencies := list, parameters := rec( galoisField:=galoisField, C:=C, D:=D, i:=number )) ;
end;

GAUSS_ExtendParameters := function(i, j, TaskListClearDown, TaskListE)
	local list, A, E, flag; # parameters for GAUSS_Extend as in subprograms.g
	
    if (j = 1) then
		list := [ TaskListClearDown[i][j] ];
		E := rec( rho := [], delta := [], nr := 0 );
    else
		list := [ TaskListClearDown[i][j], TaskListE[i][j-1] ];
		E := TaskResult( TaskListE[i][j-1] );
    fi;
    
	A := TaskResult( TaskListClearDown[i][j] ).A;
	flag := j;
	return rec( dependencies := list, parameters := rec( A:=A, E:=E, flag:=flag) );
end;

GAUSS_UpdateRowParameters := function(i, j, k, matrixC, TaskListClearDown, TaskListUpdateR, galoisField)
	local list, A, C, B, number; # parameters for GAUSS_UpdateRow as in subprograms.g

	list := [ TaskListClearDown[i][j] ];

    if (j = 1) then
		C := matrixC[i][k];
    else
		C := TaskResult( TaskListUpdateR[i][j-1][k] ).C;
        Add(list, TaskListUpdateR[i][j-1][k]);
    fi;

    if (i = 1) then
		B := [];
    else
		B := TaskResult( TaskListUpdateR[i-1][j][k] ).B;
        Add(list, TaskListUpdateR[i-1][j][k]);
    fi;

	A := TaskResult( TaskListClearDown[i][j] ).A;
	number := i;
	return  rec( dependencies:=list, parameters:=rec(galoisField:=galoisField, A:=A, C:=C, B:=B, i:=number) );
end;

GAUSS_UpdateRowTrafoParameters := function(i, j, k, TaskListClearDown, TaskListE, TaskListUpdateM, galoisField)
	local list, A, K, M, E; # parameters for GAUSS_UpdateRowTrafe as in subprograms.g

    Info(InfoGauss, 4, "Start UpdateRowTrafoParameters", i, " ", j, " ", k);
	list := [ TaskListClearDown[i][j], TaskListE[k][j] ];
    
    if (i = 1) then
        M := [];
    else
		M := TaskResult( TaskListUpdateM[i-1][j][k] ).M;
        Add(list, TaskListUpdateM[i-1][j][k]); 
    fi;

    if (j = 1) then
		K := [];
    else
		K := TaskResult( TaskListUpdateM[i][j-1][k] ).K;
        Add(list, TaskListUpdateM[i][j-1][k]);
    fi;
    
    A := TaskResult( TaskListClearDown[i][j] ).A;
    E := TaskResult( TaskListE[k][j] );

	return rec( dependencies:=list, parameters:=rec(galoisField:=galoisField, A:=A, K:=K, M:=M, E:=E, i:=i, k:=k, j:=j) );
end;
