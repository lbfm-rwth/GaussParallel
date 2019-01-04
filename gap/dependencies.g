# This file contains the help functions that manage the parallelization of the
# algorithm. In particular, the functions in this file manage the dependencies
# of the algorithm's subprograms between each other.

# This function is temporary until it is not necessary anymore to find
# a number of blocks that divides the width or height.
# When different chop sizes are possible, we will presumably find the perfect
# chop size and divide in blocks of this size.
GAUSS_calculateBlocks := function( a )
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

    if not (i = 1) then
        Add(list, TaskListClearDown[i-1][j]);
    fi;

    if not (j = 1) then
        Add(list, TaskListUpdateR[i][j-1][j]);
    fi;

    return rec( dependencies := list );
end;

GAUSS_ExtendParameters := function(i, j, TaskListClearDown, TaskListE)
    local list, A, E, flag; # parameters for GAUSS_Extend as in subprograms.g

    if (j = 1) then
        list := [ TaskListClearDown[i][j] ];
    else
        list := [ TaskListClearDown[i][j], TaskListE[i][j-1] ];
    fi;

    return rec( dependencies := list );
end;

GAUSS_UpdateRowParameters := function(i, j, k, matrixC, TaskListClearDown, TaskListUpdateR, galoisField)
    local list, A, C, B, number; # parameters for GAUSS_UpdateRow as in subprograms.g

    list := [ TaskListClearDown[i][j] ];

    if not (j = 1) then
        Add(list, TaskListUpdateR[i][j-1][k]);
    fi;

    if not (i = 1) then
        Add(list, TaskListUpdateR[i-1][j][k]);
    fi;

    return  rec( dependencies:=list );
end;

GAUSS_UpdateRowTrafoParameters := function(i, j, k, TaskListClearDown, TaskListE, TaskListUpdateM, galoisField)
    local list, A, K, M, E; # parameters for GAUSS_UpdateRowTrafe as in subprograms.g

    Info(InfoGauss, 4, "Start UpdateRowTrafoParameters", i, " ", j, " ", k);
    list := [ TaskListClearDown[i][j], TaskListE[k][j] ];

    if not (i = 1) then
        Add(list, TaskListUpdateM[i-1][j][k]);
    fi;

    if not (j = 1) then
        Add(list, TaskListUpdateM[i][j-1][k]);
    fi;

    return rec( dependencies:=list );
end;
