HandlePivots := function( j,n,p,Eold )
    local v,w,E;

    if j = 1 then       
        E := 0*[ 1 .. n ];
    fi;

           
    if not j = 1 then 
        E := ExtendPivotRows( Eold,p );
        v := MKR( Eold,E );
        w := MKw( Eold,E );
    else    
        v := E;
        E := ExtendPivotRows( E,p );
        v := MKR(v,E );
        w := E;
    fi;

    return [E,v,w];
end;

MAD := function( X,Y,Z )
    if IsEmpty(Z) then
        return [];
    fi;
    return X + Y*Z;
end;

GaussParallelTrafo := function( Inp,a,b,f ) #Chop inputmatrix Inp into (a)x(b) matrix
    local C,D,B,A,E,F,M,K,X,R, tmp,tmpR,tmpC,i,j,k,h,v,w,rank,rows,
    dummyTask, TaskListClearDown,TaskListUpdateRow,TaskListPivots,
    TaskListUpdateTrafo,TaskListClearUp,TaskListClearUpTrafo ,ncols;
    C := ChoppedMatrix( f,Inp,a,b );w := []; v :=[];R := [];A := [];B := [];D := [];E := [];F := [];M := [];K := [];X := [];
    ncols := DimensionsMat( Inp )[2];
    Inp := MutableCopyMat(C);   
    
    dummyTask := RunTask( function() return [ [],[]  ]; end  ); 
    # Initialisation of data sets
    TaskListClearDown := List(
        [1..a],
        x -> List( [1..b], x -> dummyTask )
    );
    TaskListClearUpTrafo := List(
        [1..b],
        x -> List(
            [1..a],
            x -> List( [1..b], x -> dummyTask )
        )
    );
    TaskListClearUp := List(
        [1..b],
        x -> List(
            [1..b],
            x -> List( [1..b], x -> dummyTask )
        )
    );

    TaskListPivots := List(
        [1..a],
        x -> List( [1..b], x -> dummyTask )
    );
    TaskListUpdateRow := List(
        [1..a],
        x -> List(
            [1..b],
            x -> List( [1..b], x -> dummyTask )
        )
    );
    TaskListUpdateTrafo := List(
        [1..a],
        x -> List(
            [1..b],
            x -> List( [1..b], x -> dummyTask )
        )
    );

    for i in [ 1 .. a ] do
        A[i] := [];
        E[i] := [];
	K[i] := [];
        v[i] := [];
        w[i] := [];
        for k in [ 1 .. b ] do
            A[i][k] := [];
            E[i][k] := [];
            v[i][k] := [];
            w[i][k] := [];
        od;
        for h in [ 1 .. a ] do
            K[i][h] := [];
        od;
    od;
    for i in [ 1 .. b ] do
        M[i] := [];
        F[i] := [];
        for k in [ 1 .. a ] do
            M[i][k] := [];
            F[i][k] := [];
        od;
    od;
    for k in [ 1 .. b ] do
        D[k] := rec(remnant := [], pivots := []);
        B[k] := [];
        R[k] := [];
        X[k] := [];
        for j in [ 1 .. b ] do
            B[k][j] := [];
            R[k][j] := [];
            X[k][j] := [];
        od;
    od;
  
    # Step1: ClearDown
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            if i=1 and j=1 then
                TaskListClearDown[1][1] :=
                RunTask( ClearDown,f,Inp[1][1],[],[] );
            elif i=1 and j>1 then
                TaskListClearDown[i][j] := ScheduleTask(
                    TaskListUpdateRow[i][j-1][j],
                    ClearDown, 
                    f, 
                    TaskResult( TaskListUpdateRow[i][j-1][j] )[1],
                    [],
                    []
                );
            elif i>1 and j=1 then
                TaskListClearDown[i][j] := ScheduleTask(
                    TaskListClearDown[i-1][j],
                    ClearDown, 
                    f, 
                    Inp[i][j],
                    TaskResult( TaskListClearDown[i-1][j] )[2],
                    TaskResult( TaskListClearDown[i-1][j] )[1]
                );
            else
                TaskListClearDown[i][j] := ScheduleTask(
                    [ TaskListClearDown[i-1][j],
                      TaskListUpdateRow[i][j-1][j] ],
                    ClearDown, 
                    f, 
                    TaskResult( TaskListUpdateRow[i][j-1][j] )[1],
                    TaskResult( TaskListClearDown[i-1][j] )[2],
                    TaskResult( TaskListClearDown[i-1][j] )[1]
                );
            fi;

            if j=1 then
            TaskListPivots[i][j] := ScheduleTask(
                [ TaskListClearDown[i][j] ],
                HandlePivots,
                j,
                DimensionsMat(Inp[i][j])[1],
                TaskResult(TaskListClearDown[i][j])[3][5],
                []
            ); 
            else
            TaskListPivots[i][j] := ScheduleTask(
                [ TaskListPivots[i][j-1],
                TaskListClearDown[i][j] ],
                HandlePivots,
                j,
                DimensionsMat(Inp[i][j])[1],
                TaskResult(TaskListClearDown[i][j])[3][5],
                TaskResult(TaskListPivots[i][j-1])[1]
            ); 
            fi;
        
        

            for k in [ j+1 .. b ] do

               if i=1 and j=1 then
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        [ TaskListClearDown[i][j] ],
                        UpdateRow,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        Inp[i][k],
                        []
                    );
                elif i=1 and j>1 then
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateRow[i][j-1][k] ],
                        UpdateRow,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        TaskResult( TaskListUpdateRow[i][j-1][k] )[1],
                        []
                    );
                elif i>1 and j=1 then
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateRow[i-1][j][k] ],
                        UpdateRow,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        Inp[i][k],
                        TaskResult( TaskListUpdateRow[i-1][j][k] )[2]
                    );
                else
                    TaskListUpdateRow[i][j][k] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateRow[i][j-1][k],
                          TaskListUpdateRow[i-1][j][k] ],
                        UpdateRow,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        TaskResult( TaskListUpdateRow[i][j-1][k] )[1],
                        TaskResult( TaskListUpdateRow[i-1][j][k] )[2]
                    );
                fi;
            od;
        
            for h in [ 1 .. i-1 ] do 
               if j=1 then
                    TaskListUpdateTrafo[i][j][h] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateTrafo[i-1][j][h], 
                          TaskListPivots[h][j] ],
                        UpdateTrafo,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        [],
                        TaskResult( TaskListUpdateTrafo[i-1][j][h] )[2],
                        TaskResult( TaskListPivots[h][j] )[2],
                        0,
                        TaskResult( TaskListPivots[h][j] )[3]
                    );
                else
                    TaskListUpdateTrafo[i][j][h] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateTrafo[i][j-1][h],
                          TaskListUpdateTrafo[i-1][j][h], 
                          TaskListPivots[h][j] ],
                        UpdateTrafo,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        TaskResult( TaskListUpdateTrafo[i][j-1][h] )[1],
                        TaskResult( TaskListUpdateTrafo[i-1][j][h] )[2],
                        TaskResult( TaskListPivots[h][j] )[2],
                        0,
                        TaskResult( TaskListPivots[h][j] )[3]
                    );
                fi;
            od;

            if i=1 and j=1 then
                    TaskListUpdateTrafo[i][j][i] := ScheduleTask(
                        [ TaskListClearDown[i][j], 
                          TaskListPivots[i][j] ],
                        UpdateTrafo,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        [],
                        [],
                        TaskResult( TaskListPivots[i][j] )[2],
                        1,
                        TaskResult( TaskListPivots[i][j] )[3]
                    );
                         
            elif i=1 and j>1 then
                    TaskListUpdateTrafo[i][j][i] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateTrafo[i][j-1][i], 
                          TaskListPivots[i][j] ],
                        UpdateTrafo,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        TaskResult( TaskListUpdateTrafo[i][j-1][i] )[1],
                        [],
                        TaskResult( TaskListPivots[i][j] )[2],
                        1,
                        TaskResult( TaskListPivots[i][j] )[3]
                    );

            elif i>1 and j=1 then
                    TaskListUpdateTrafo[i][j][i] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateTrafo[i-1][j][i], 
                          TaskListPivots[i][j] ],
                        UpdateTrafo,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        [],
                        TaskResult( TaskListUpdateTrafo[i-1][j][i] )[2],
                        TaskResult( TaskListPivots[i][j] )[2],
                        1,
                        TaskResult( TaskListPivots[i][j] )[3]
                    );

            else
                    TaskListUpdateTrafo[i][j][i] := ScheduleTask(
                        [ TaskListClearDown[i][j],
                          TaskListUpdateTrafo[i][j-1][i],
                          TaskListUpdateTrafo[i-1][j][i], 
                          TaskListPivots[i][j] ],
                        UpdateTrafo,
                        f,
                        TaskResult( TaskListClearDown[i][j] )[3],
                        TaskResult( TaskListUpdateTrafo[i][j-1][i] )[1],
                        TaskResult( TaskListUpdateTrafo[i-1][j][i] )[2],
                        TaskResult( TaskListPivots[i][j] )[2],
                        1,
                        TaskResult( TaskListPivots[i][j] )[3]
                    );

            fi;
        od;
       
        # F
        for j in [ 1 .. b ] do
            F[j][i] := ScheduleTask(
            TaskListPivots[i][b],
            MKR,
            TaskResult(TaskListPivots[i][j])[1],
            TaskResult(TaskListPivots[i][b])[1]
            );
        od;
    od;
    Print( "tasks succesfully scheduled\n" );
    WaitTask( Concatenation( TaskListClearDown ) );
    WaitTask( Concatenation( F ) );
    WaitTask( Concatenation( TaskListPivots ) );
    WaitTask( Concatenation( List(TaskListUpdateRow,Concatenation) ) );
    WaitTask( Concatenation( List(TaskListUpdateTrafo,Concatenation) ) );

    #Error("CHECK ALL TASKRESULTS");
    # convert to old output format for simplicity first - to be chenged

    for i in [ 1 .. a ] do
        for k in [ 1 .. b ] do
            A[i][k] := TaskResult(TaskListClearDown[i][k])[3];
            E[i][k] := TaskResult(TaskListPivots[i][k])[1];
            v[i][k] := TaskResult(TaskListPivots[i][k])[2];
            w[i][k] := TaskResult(TaskListPivots[i][k])[3];
            M[k][i] := TaskResult(TaskListUpdateTrafo[a][k][i])[2];
        od;
    od;
    for k in [ 1 .. b ] do
        D[k].pivots := TaskResult(TaskListClearDown[a][k])[2];
        D[k].remnant := TaskResult(TaskListClearDown[a][k])[1];
        for j in [ 1 .. b ] do
            B[k][j] := TaskResult( TaskListUpdateRow[a][k][j] )[2];
        od;
    od;

    # Step2: Riffle missing collumns into the M_jh's
    for j in [ 1 .. b ] do
        for h in [ 1 .. a ] do
            M[j][h] := RowLenghten( f,M[j][h],TaskResult(F[j][h])  );
        od;
    od;

    Print("Upwards Cleaning...");

    # Step3: Upwards Cleaning
    for i in [ 1 .. b ] do #we use i for b-k+1 from now on
       k := b-i+1;
       R[k][k] := D[k].remnant;
    od;

    for i in [ 1 .. b ] do
        k := b-i+1;

        for j in [ 1 .. k-1 ] do
            tmp := CEX( f,BitstringToCharFct(D[k].pivots),B[j][k] );
            # ConvertToMatrixRepNC( tmp[1],f );
            # ConvertToMatrixRepNC( tmp[2],f );
         
            X[j][k] := tmp[1];
            R[j][k] := tmp[2];
        od;
    od; 
    
    for i in [ 1 .. b ] do
        k := b-i+1;
        for j in [ 1 .. k-1 ] do
            for h in [ k .. b ] do
                if h-k=0 then
                    TaskListClearUp[j][h][h-k+1] := ScheduleTask(
                        [],
                        MAD,
                        R[j][h],
                        X[j][k],
                        R[k][h]
                    );
                else
                    TaskListClearUp[j][h][h-k+1] := ScheduleTask(
                        [TaskListClearUp[k][h][h-k],
                        TaskListClearUp[j][h][h-k]],
                        MAD,
                        TaskResult(TaskListClearUp[j][h][h-k]),
                        X[j][k],
                        TaskResult(TaskListClearUp[k][h][h-k])
                    );
                fi;
            od;
            for h in [ 1 .. a ] do
                if i=1 then
                    
                    TaskListClearUpTrafo[j][h][i] := ScheduleTask(
                        [],
                        MAD,
                        M[j][h],
                        X[j][k],
                        M[k][h]
                    ); 
                else
                    TaskListClearUpTrafo[j][h][i] := ScheduleTask(
                        [TaskListClearUpTrafo[k][h][i-1],
                        TaskListClearUpTrafo[j][h][i-1]],
                        MAD,
                        TaskResult(TaskListClearUpTrafo[j][h][i-1]),
                        X[j][k],
                        TaskResult(TaskListClearUpTrafo[k][h][i-1])
                    );
                fi;
            od;
        od;
    od;

    Print( "tasks succesfully scheduled\n" );

    WaitTask( Concatenation( List(TaskListClearUp,Concatenation) ) );
    WaitTask( Concatenation( List(TaskListClearUpTrafo,Concatenation) ) );

    for i in [ 1 .. b ] do
        k := b-i+1;
        for j in [ 1 .. k-1 ] do
            for h in [ k .. b ] do
                #if not h-k=0 then
                    R[j][h] := TaskResult(TaskListClearUp[j][h][h-k+1]);
                #fi;
            od;
            for h in [ 1 .. a ] do
                #if not i=1 then
                    M[j][h] := TaskResult(TaskListClearUpTrafo[j][h][i]);
                #fi;
            od;
        od;
    od;


       # for j in [ 1 .. k-1 ] do
       #    tmp := CEX( f,BitstringToCharFct(D[k].pivots),B[j][k] );
       #   # ConvertToMatrixRepNC( tmp[1],f );
       #   # ConvertToMatrixRepNC( tmp[2],f );
       #    
       #    X[j][k] := tmp[1];
       #    R[j][k] := tmp[2];
       # od; 

        #   for h in [ k .. b ] do
        #      if not IsEmpty(R[k][h]) then
        #       
        #         ConvertToMatrixRepNC( R[j][h],f );
        #         ConvertToMatrixRepNC( R[k][h],f );
        #         R[j][h] := R[j][h] + X[j][k]*R[k][h];
        #      fi;
        #   od;
        #   for h in [ 1 .. a ] do
        #      if not IsEmpty(M[k][h]) then
        #         ConvertToMatrixRepNC( M[j][h],f );
        #         ConvertToMatrixRepNC( M[k][h],f );
        #         M[j][h] := M[j][h] + X[j][k]*M[k][h];
        #      fi;
        #   od;
        #od;


    ## Write output

    # Beginning with row-select bitstring now named v
    v := [];
    rank := 0;
    for i in [ 1 .. a ] do
        v := Concatenation( v,E[i][b] );
        rank := rank + Sum( E[i][b] );
    od;

    # Glueing the blocks of M
    B := NullMat( rank,Length(v),f );
    rows := [];

    for j in [ 1 .. b ] do
        rows[j] := 0;
        for i in [ 1 .. a ] do
            if  not IsEmpty(M[j][i]) then
                rows[j] := DimensionsMat(M[j][i])[1];
                break;
            fi;
        od;
    od;

    tmpR := 1;
    for j in [ 1 .. b ] do
        if rows[j]=0 then
            continue;
        fi;
        tmpC := 1;
        for i in [ 1 .. a ] do
            if IsEmpty(M[j][i]) then
                M[j][i] := NullMat( rows[j],Length(E[i][b]),f );
            else
                M[j][i] := TransposedMat( RRF( f,NullMat( Length(E[i][b])-DimensionsMat(M[j][i])[2],DimensionsMat(M[j][i])[1],f ),TransposedMat(M[j][i]),E[i][b] ) );
            fi;

            B{[tmpR .. tmpR + DimensionsMat(M[j][i])[1]-1 ]}{[tmpC .. tmpC + DimensionsMat(M[j][i])[2]-1 ]}
            := M[j][i];
            tmpC := tmpC + DimensionsMat(M[j][i])[2];
        od;
        tmpR := tmpR + rows[j];
    od;

    # Glueing the blocks of R
    #Error("BREAKPOINT - Checking remnant"); 
    C := NullMat( rank,ncols-rank,f );
    rows := [];
    w := [];

    for i in [ 1 .. b ] do
         rows[i] := 0;
         if IsEmpty(D[i].pivots) then
             w := Concatenation( w,0*[1..DimensionsMat(Inp[1][i])[2]] );
         else
             w := Concatenation( w,D[i].pivots );
         fi;
         for j in [ 1 .. b ] do
             if  not IsEmpty(R[i][j]) then
                 rows[i] := DimensionsMat(R[i][j])[1];
                 break;
             fi;
         od;
     od;
 
     tmpR := 1;
     for i in [ 1 .. b ] do
         if rows[i]=0 then
             continue;
         fi;
         tmpC := 1;
         for j in [ 1 .. b ] do
             if IsEmpty(R[i][j]) then
                 if not IsEmpty(D[j].pivots) then
                     tmpC := tmpC + Sum( 1 - D[j].pivots );
                 elif  not IsEmpty(R[1][j]) then
                     tmpC := tmpC + DimensionsMat(R[1][j])[2];
                 fi;
                 continue;            
             fi;
 
             C{[tmpR .. tmpR + DimensionsMat(R[i][j])[1]-1 ]}{[tmpC .. tmpC + DimensionsMat(R[i][j])[2]-1 ]}
             := R[i][j];
             tmpC := tmpC + DimensionsMat(R[i][j])[2];
         od;
         tmpR := tmpR + rows[i];
     od;    

     # SLOW - only for testing 
     C := TransposedMat( RRF( f,TransposedMat(C), -IdentityMat( rank,f ),w  ) );

     return [v,w,C,B,D];
end;
