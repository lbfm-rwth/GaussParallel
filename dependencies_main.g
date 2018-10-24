# This file contains some help functions making the parallel parts of the code 
# in main.g more readable. (Dependency of main.g.)

GAUSS_ClearUp := function( R,X,R_ )
    if IsEmpty(R_) or IsEmpty(X) then return R; fi;
    if IsEmpty(R) then
        return X*R_;
    else
        return R + X*R_;
    fi;
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
	return [ list, galoisField, C, D, number ];
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
	return [ list, A, E, flag ];
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
	return [ list, galoisField, A, C, B, number ];
end;

GAUSS_UpdateRowTrafoParameters := function(i, j, k, TaskListClearDown, TaskListE, TaskListUpdateM, galoisField)
	local list, A, K, M, E; # parameters for GAUSS_UpdateRowTrafe as in subprograms.g

    Info(InfoGauss, 3, "Start UpdateRowTrafoParameters", i, " ", j, " ", k);
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

	return [ list, galoisField, A, K, M, E, i, k, j];
end;
