###############################
# function ClearDown
# Input:
#   f - Field
#   i,j - where to find the following:
#       H - Input matrix
#       t - Bitlist of pivot columns positions
#       R - Residue of above block after echolonization
#
# Output:
#   [ R, t, T ]
#   R - ??
#   t - ??
#   T - ??
###############################
ClearDown := function( f, H, t, R )
    local tmp, Chi, ct, HH, tt, RR, ttt, RRR, i, RRn, A, AA, T, M, K, E, s, u;

    H := ShallowCopy( H );
    t := ShallowCopy( t );
    R := ShallowCopy( R );

   # if not Length(t) = DimensionsMat(H)[2] then
   #     Error( "Length of bitlist t does not match dimensions of H!" );
   # fi;

    #### INITIALIZATION ####
    ## Residue R was empty. Thus matrix above has full column-rank.
    if not IsEmpty(t) and IsEmpty( R ) then
        A := H;
        M := [];
        E := [];
        K := [];
        s := [];
        RR := [];
        tt := [];
        ttt := t;
        u := 0 * t;
        T:=[A, M, E, K, s, u];
        return [RR, ttt, T];
    fi;

    ## First step (i=1,  i block row) or all matrices above have rank 0
    if IsEmpty( t ) then
        # A is empty iff this case happens, Q is empty then aswell
        A := [];
        HH := H;
    else
        # Column Extraction
        tmp := CEX( f, BitstringToCharFct(t), H );
        A := tmp[1]; AA := tmp[2];
        # Reduce H to (0|HH)
        # Mult Add
        HH := AA + A*R;
    fi;
 
 
    #### END INITIALIZATION ####

    # Echelonization
    tmp := ECH( f, HH );
    M:=tmp[1];K:=tmp[2];RR:=tmp[3];s:=tmp[4];tt:=tmp[5];
    #Error( "Break Point - echel" );

    # TODO complement then extend?
    if IsEmpty(t) then Chi := tt;
    else
     Chi := 0*[1..DimensionsMat(H)[2]];ct:=1;
     if not tt=[] then
        for i in [1..Length(Chi)] do
            if t[i]=0 then
                if tt[ct]=1 then Chi[i] := 1; fi;
                ct:= ct+1;
            fi;
        od;
     fi;
    fi;

    ## Are we running into any special cases, where return values
    ## of the above echelonization are empty?
    # The case K empty is handled by UpdateRow.
    #
    # The following are equivalent:
    # s is empty
    # tt is empty
    # M is empty
    #
    # This only happens when A*R + AA = 0
    if ForAny( [ M, s, tt ], IsEmpty ) then
        M := [];
        E := [];
        K := [];
        s := [];
        RR := R; # Residue does not change
        ttt := t;
        u := 0 * t;
        T:=[A, M, E, K, s, u];
        return [RR, ttt, T];
    fi;
    # If RR is empty, but tt is not, then the bitstring tt, representing
    # the positions of the new pivot columns, is AllOne.
    # In this case, there is nothing to be done here.

    #Error( "Break Point - before CEX new residue" );
    tmp := CEX( f, BitstringToCharFct(tt), R );
    E:=tmp[1];RRn:=tmp[2];
    ## Update the residue and the pivot column bitstring
    tmp := PVC( BitstringToCharFct(t), BitstringToCharFct(Chi) );
    ttt:=CharFctToBitstring(DimensionsMat(H)[2], tmp[1]); u:=tmp[2];
    # Error( "Break Point - after CEX new residue" );
    
    T:=[A, M, E, K, s, u];

    ## Did column extraction return empty values?
    if IsEmpty(E) then ## if the above was all zero but we got new pivots in the current iteration
        return [RR, ttt, T];
    fi;

    ## RRn is empty, iff. the new pivot columns completely
    ## annihilate the old residue.
    if IsEmpty(RRn) then
        RR := [];
    else
        RRR:=RRn+E*RR;
        RR := RRF( RRR, RR, u );
    fi;

    return [RR, ttt, T];
end;

UpdateRow := function( f, T, H, Bjk )
 local A, E, M, K, s, u,  tmp, Z, V, X, W, S, B;
 T := ShallowCopy(T);
 H := ShallowCopy(H);
 B := ShallowCopy(Bjk);
 A:=T[1];M:=T[2];E:=T[3];K:=T[4];s:=T[5];u:=T[6];
 
 ###
 # If A is empty, there are no rowoperations form above to consider
 ###
 if IsEmpty(A) then
  Z := H;
 else 
  Z := A*B+H;
 fi;

 tmp := REX( f, BitstringToCharFct(s), Z );
 V:=tmp[1];W:=tmp[2];
 ###
 # If V is empty, then there where no operations exept from A
 # in this case there is nothing more to update
 ###
 if IsEmpty(V) then
  return [Z,B]; 
 else 
  X:=M*V;
 fi;

 S:= E*X+B;
 B:=RRF( S, X, u );
 
 ###
 # if K is empty, then s is the all-one-bitstring and there are no non-pivot rows
 # which would change according to K. So K should be empty and there is nothing more to update
 ###
 if not IsEmpty(K) then
  # s is neither empty nor all-one at this point
  H := K*V+W;
 else
  H := W;
 fi;
 
 return [H, B];
end;

Step1 := function( A, n )
    local f, C, B, Rj, tj,
        dummyTask, TaskListClearDown, TaskListUpdateRow,
        i, j, k;
    ## Chop A into an nxn matrix
    f := DefaultFieldOfMatrix( A );
    C := ChopMatrix( n, A );
    # FIXME UNUSED CODE
    ## Initialize B as an nxn list pointing to empty lists
    B := List( [1..n], x -> List( [1..n], x -> [] ) );
    Rj := [];
    tj := [];
    ########## Keep track of Tasks ##########
    dummyTask := RunTask( function() return 0; end );
    ## TODO Do we need to erase the TaskResults manually?
    ## nxn many tasks
    TaskListClearDown := List(
        [1..n],
        x -> List( [1..n], x -> dummyTask )
    );
    ## nxnxn many tasks
    TaskListUpdateRow := List(
        [1..n],
        x -> List(
            [1..n],
            x -> List( [1..n], x -> dummyTask )
        )
    );
    for i in [ 1 .. n ] do
        for j in [ 1 .. n ] do
            ########## Schedule ClearDown Tasks ##########
            ## first row, first column: start computation
            if i = 1 and j = 1 then
                TaskListClearDown[i][j] := RunTask(
                    function()
                        local H, t, R;
                        Print("Starting computation!");
                        Error( "Break Point - First Task!" );
                        H := C[i][j];
                        t := [];
                        R := [];
                        return ClearDown( f, H, t, R );
                    end
                );
            ## first row: wait for left-side `UpdateRow`s
            elif i = 1 and j > 1 then
                TaskListClearDown[i][j] := ScheduleTask(
                    ## Condition
                    TaskListUpdateRow[i][j-1][j],
                    ## Function Call
                    function()
                        local H, t, R;
                        H := TaskResult( TaskListUpdateRow[i][j-1][j] )[1];
                        t := [];
                        R := [];
                        return ClearDown( f, H, t, R );
                    end
                );
            ## first column: wait for upper `ClearDown`s
            elif i > 1 and j = 1 then
                TaskListClearDown[i][j] := ScheduleTask(
                    ## Condition
                    TaskListClearDown[i-1][j],
                    ## Function Call
                    function()
                        local H, t, R;
                        H := C[i][j];
                        t := TaskResult( TaskListClearDown[i-1][j] )[2];
                        R := TaskResult( TaskListClearDown[i-1][j] )[1];
                        return ClearDown( f, H, t, R );
                    end
                );
            else ## i > 1, j > 1: wait for "everything"
                TaskListClearDown[i][j] := ScheduleTask(
                    ## Condition: List of tasks to wait on
                    [ TaskListClearDown[i-1][j],
                      TaskListUpdateRow[i][j-1][j] ],
                    ## Function Call
                    function()
                        local H, t, R;
                        H := TaskResult( TaskListUpdateRow[i][j-1][j] )[1];
                        t := TaskResult( TaskListClearDown[i-1][j] )[2];
                        R := TaskResult( TaskListClearDown[i-1][j] )[1];
                        return ClearDown( f, H, t, R );
                    end
                );
            fi;
            ########## Schedule UpdateRow Tasks ##########
            for k in [ j+1 .. n ] do
                ## first row: since j = 2 no previous UpdateRow was spawned
                if i = 1 and j = 1 then
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        ## Condition: List of tasks to wait on
                        [ TaskListClearDown[i][j] ],
                        ## Function Call
                        function()
                            local T, H, B;
                            T := TaskResult( TaskListClearDown[i][j] )[3];
                            H := C[i][k];
                            B := [];
                            return UpdateRow( f, T, H, B );
                        end
                    );
                ## first row: wait
                elif i = 1 and j > 1 then
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        ## Condition: List of tasks to wait on
                        [ TaskListClearDown[i][j],
                          TaskListUpdateRow[i][j-1][k] ],
                        ## Function Call
                        function()
                            local T, H, B;
                            T := TaskResult( TaskListClearDown[i][j] )[3];
                            H := TaskResult( TaskListUpdateRow[i][j-1][k] )[1];
                            B := [];
                            return UpdateRow( f, T, H, B );
                        end
                    );
                elif i > 1 and j = 1 then
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        ## Condition: List of tasks to wait on
                        [ TaskListClearDown[i][j],
                          TaskListUpdateRow[i-1][j][k] ],
                        ## Function Call
                        function()
                            local T, H, B;
                            T := TaskResult( TaskListClearDown[i][j] )[3];
                            H := C[i][k];
                            B := TaskResult( TaskListUpdateRow[i-1][j][k] )[2];
                            return UpdateRow( f, T, H, B );
                        end
                    );
                else ## i > 1 and j > 1
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        ## Condition: List of tasks to wait on
                        [ TaskListClearDown[i][j],
                          TaskListUpdateRow[i][j-1][k],
                          TaskListUpdateRow[i-1][j][k] ],
                        ## Function Call
                        function()
                            local T, H, B;
                            T := TaskResult( TaskListClearDown[i][j] )[3];
                            H := TaskResult( TaskListUpdateRow[i][j-1][k] )[1];
                            B := TaskResult( TaskListUpdateRow[i-1][j][k] )[2];
                            return UpdateRow( f, T, H, B );
                        end
                    );
                fi;
            od;
        od;
    od;
    ## TODO V is that so? V
    ## This is implicitly waiting on all UpdateRow calls
    WaitTask( Concatenation( TaskListClearDown ) );
    WaitTask( Concatenation( List( TaskListUpdateRow, Concatenation ) ) );
    tj := List( [ 1..n ], j -> TaskResult( TaskListClearDown[n][j] )[2] );
    Rj := List( [ 1..n ], j -> TaskResult( TaskListClearDown[n][j] )[1] );
    Error( "Break Point - END OF STEP1" );
    return [ C, B, Rj, tj ];
end;

_DEPRECATED_Step1 := function( A,n )
 local C, cur,f, nrr, tmp,   i, j, k, B, T, Rj, H, tj, V, W;
 f := DefaultFieldOfMatrix( A );
 C := ChopMatrix( n, A );
 nrr := DimensionsMat(A)[1];
 B := []; #Init B as nxn; 
 for i in [1..n] do
  B[i] :=[];  
  for j in [1..n] do
   B[i][j]:=[];
  od;
 od;
 Rj := [];
 tj := [];

 for i in [1..n] do
  for j in [1..n] do
   cur := 1;
   H := C[i][j];

   if i = 1 then
    tj[j] := [];
    Rj[j] := [];
   fi;
  
   if IsEmpty(H) then continue; fi;
   tmp := ClearDown( f, H, tj[j], Rj[j] );
   Rj[j] := tmp[1];
   tj[j] := tmp[2];
   T := tmp[3];
   
   for k in [j+1..n] do
    if i = 1 then
      B[j][k]:=[];
    fi; 
    
    tmp := UpdateRow( f, T, C[i][k], B[j][k] );
    C[i][k] := tmp[1];
    B[j][k] := tmp[2];
   od;
  od;
 od;

 return [ C, B, Rj, tj];
end;
