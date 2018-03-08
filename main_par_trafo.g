ChiefParallelClearDown := function( galoisField,mat,a,b )
    ## inputs: a finite field, a matrix, natural nr. a,b to treat mat as axb matirx
    local   TaskListClearDown,
            TaskListUpdateR,
            TaskListUpdateM,
            TaskListE,
            dummyTask,
            nrows,
            ncols,
            rows,
            rank,
            tmpC,
            tmpR,
            tmp,
            bs,
            A,
            B,
            C,
            D,
            E,
            K,
            M,
            R,
            X,
            v,
            w,
            i,
            j,
            k,
            k_,
            l,
            h;

    ##Preparation: Init and chopping the matrix mat
    C := ChopMatrix( galoisField,mat,a,b );
    ncols := DimensionsMat( mat )[2];
    dummyTask := RunTask( function() return []; end );
    TaskListClearDown := List(
        [ 1 .. a ],
        x -> List( [ 1 .. b ], x -> dummyTask )
    );
    TaskListE := List(
        [ 1 .. a ],
        x -> List( [ 1 .. b ], x -> dummyTask )
    );
    TaskListUpdateR := List(
        [ 1 .. b ],
        x -> List( 
            [ 1 .. b ], 
            x -> List( [ 1 .. b ], x -> RunTask( function() return rec( C := [],B:=[] ); end ) )
        )
    );
    TaskListUpdateM := List(
        [ 1 .. b ],
        x -> List( 
            [ 1 .. a ], 
            x -> List( [ 1 .. b ], x -> RunTask( function() return rec( M := [],K:=[] ); end ) )
        )
    );

    A := [];
    B := [];
    D := [];
    E := [];
    K := [];
    M := [];
    R := [];
    X := [];
    nrows := [];
    for i in [ 1 .. a ] do
        A[i] := [];
        E[i] := [];
        K[i] := [];
        for k in [ 1 .. b ] do
            A[i][k] := rec( A:=[],M:=[],K:=[],rho:=[],E:=[],lambda:=[] );
            E[i][k] := rec( rho:=[],delta:=[],nr:=0 );
        od;
        for h in [ 1 .. a ] do
            K[i][h] := [];
        od;
    od;
    for i in [ 1 .. b ] do
        M[i] := [];
        nrows[i] := DimensionsMat(C[1][i])[2];
        for k in [ 1 .. a ] do
            M[i][k] := [];
        od;
    od;
    for k in [ 1 .. b ] do
        D[k] := rec( vectors:=[],bitstring:=[] );
        B[k] := [];
        R[k] := [];
        X[k] := [];
        for j in [ 1 .. b ] do
            B[k][j] := [];
            R[k][j] := [];
            X[k][j] := [];;
        od;
    od;
    ###############################
    ###############################

    ## Step 1 ##
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            if i=1 and j=1 then
                TaskListClearDown[1][1] := RunTask(
                    ClearDown,
                    galoisField,
                    C[1][1],
                    rec( vectors:=[],bitstring:=[] ),
                    1
                );
                TaskListE[i][j] := ScheduleTask(
                    [
                        TaskListClearDown[i][j]
                    ],
                    Extend,
                    TaskResult( TaskListClearDown[i][j] ).A,
                    rec( rho:=[],delta:=[],nr:=0 ),
                    j
                );
            elif i=1 and j>1 then
                TaskListClearDown[1][j] := ScheduleTask(
                    [
                        TaskListUpdateR[1][j-1][j]
                    ],
                    ClearDown,
                    galoisField,
                    TaskResult( TaskListUpdateR[1][j-1][j] ).C,
                    rec( vectors:=[],bitstring:=[] ),
                    1
                );
                TaskListE[i][j] := ScheduleTask(
                    [
                        TaskListClearDown[i][j],
                        TaskListE[i][j-1],

                    ],
                    Extend,
                    TaskResult( TaskListClearDown[i][j] ).A,
                    TaskResult( TaskListE[i][j-1] ),
                    j
                );
            elif i>1 and j=1 then
                TaskListClearDown[i][1] := ScheduleTask(
                    [
                        TaskListClearDown[i-1][1]
                    ],
                    ClearDown,
                    galoisField,
                    C[i][1],
                    TaskResult( TaskListClearDown[i-1][1] ).D,
                    i
                );
                TaskListE[i][j] := ScheduleTask(
                    [
                        TaskListClearDown[i][j]
                    ],
                    Extend,
                    TaskResult( TaskListClearDown[i][j] ).A,
                    rec( rho:=[],delta:=[],nr:=0 ),
                    j
                );
            else
                TaskListClearDown[i][j] := ScheduleTask(
                    [
                        TaskListClearDown[i-1][j],
                        TaskListUpdateR[i][j-1][j]
                    ],
                    ClearDown,
                    galoisField,
                    TaskResult( TaskListUpdateR[i][j-1][j] ).C,
                    TaskResult( TaskListClearDown[i-1][j] ).D,
                    i
                );
                TaskListE[i][j] := ScheduleTask(
                    [
                        TaskListClearDown[i][j],
                        TaskListE[i][j-1]
                    ],
                    Extend,
                    TaskResult( TaskListClearDown[i][j] ).A,
                    TaskResult( TaskListE[i][j-1] ),
                    j
                );
            fi;
           
            #tmp := ClearDown( galoisField,C[i][j],D[j],i );
            #D[j] := tmp.D;
            #A[i][j] := tmp.A;
            #if j=1 then
            #    bs := rec( rho:=[],delta:=[],nr:=0 );
            #else
            #    bs := E[i][j-1];
            #fi;
            #E[i][j] := Extend( A[i][j],bs,j );

            for k in [ j+1 .. b ] do
                if i=1 and j=1 then
                    TaskListUpdateR[i][j][k] := ScheduleTask(
                        [
                            TaskListClearDown[i][j]
                        ],
                        UpdateRow,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        C[i][k],
                        [],
                        i
                    ); 
                elif i=1 and j>1 then
                    TaskListUpdateR[i][j][k] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListUpdateR[i][j-1][k]
                        ],
                        UpdateRow,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        TaskResult( TaskListUpdateR[i][j-1][k] ).C,
                        [],
                        i
                    ); 

                elif i>1 and j=1 then
                    TaskListUpdateR[i][j][k] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListUpdateR[i-1][j][k]
                        ],
                        UpdateRow,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        C[i][k],
                        TaskResult( TaskListUpdateR[i-1][j][k] ).B,
                        i
                    );
                else
                    TaskListUpdateR[i][j][k] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListUpdateR[i-1][j][k],
                            TaskListUpdateR[i][j-1][k]
                        ],
                        UpdateRow,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        TaskResult( TaskListUpdateR[i][j-1][k] ).C,
                        TaskResult( TaskListUpdateR[i-1][j][k] ).B,
                        i
                    );
                fi;
                
                
                
                #tmp := UpdateRow(  galoisField,A[i][j],C[i][k],
                #                    B[j][k],i );
                #C[i][k] := tmp.C;
                #B[j][k] := tmp.B;
            od;
            for h in [ 1 .. i ] do
                
                if i=1 and j=1 then
                    TaskListUpdateM[i][j][h] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListE[h][j]
                        ],
                        UpdateRowTrafo,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        [],
                        [],
                        TaskResult( TaskListE[h][j] ),
                        i,
                        h,
                        j
                    ); 
                elif i=1 and j>1 then
                    TaskListUpdateM[i][j][h] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListE[h][j],
                            TaskListUpdateM[i][j-1][h]
                        ],
                        UpdateRowTrafo,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        TaskResult( TaskListUpdateM[i][j-1][h] ).K,
                        [],
                        TaskResult( TaskListE[h][j] ),
                        i,
                        h,
                        j
                    ); 

                elif i>1 and j=1 then
                    TaskListUpdateM[i][j][h] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListE[h][j],
                            TaskListUpdateM[i-1][j][h]
                        ],
                        UpdateRowTrafo,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        [],
                        TaskResult( TaskListUpdateM[i-1][j][h] ).M,
                        TaskResult( TaskListE[h][j] ),
                        i,
                        h,
                        j
                    );
                else
                    TaskListUpdateM[i][j][h] := ScheduleTask(
                        [
                            TaskListClearDown[i][j],
                            TaskListE[h][j],
                            TaskListUpdateM[i][j-1][h],
                            TaskListUpdateM[i-1][j][h]
                        ],
                        UpdateRowTrafo,
                        galoisField,
                        TaskResult( TaskListClearDown[i][j] ).A,
                        TaskResult( TaskListUpdateM[i][j-1][h] ).K,
                        TaskResult( TaskListUpdateM[i-1][j][h] ).M,
                        TaskResult( TaskListE[h][j] ),
                        i,
                        h,
                        j
                    );
                fi;
                #tmp := UpdateRowTrafo(  galoisField,A[i][j],K[i][h],
                #                   M[j][h],E[h][j],i,h,j );
                #K[i][h] := tmp.K;
                #M[j][h] := tmp.M;
            od;
        od;    
    od;

    WaitTask( Concatenation( TaskListClearDown ) );
    WaitTask( Concatenation( TaskListE ) );
    WaitTask( Concatenation( List( TaskListUpdateR,Concatenation ) ) );
    WaitTask( Concatenation( List( TaskListUpdateM,Concatenation ) ) );
    
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            E[i][j] := TaskResult( TaskListE[i][j] );
            A[i][j] := TaskResult( TaskListClearDown[i][j] ).A;
            D[j] := TaskResult( TaskListClearDown[i][j] ).D;
            for k in [ j+1 .. b ] do
                C[i][k] := TaskResult( TaskListUpdateR[i][j][k] ).C;
                B[j][k] := TaskResult( TaskListUpdateR[i][j][k] ).B;
            od;
            for h in [ 1 .. i ] do
                K[i][h] := TaskResult( TaskListUpdateM[i][j][h] ).K;
                M[j][h] := TaskResult( TaskListUpdateM[i][j][h] ).M;
            od;
        od;
    od;

    ## Step 2 ##
    for j in [ 1 .. b ] do
        for h in [ 1 .. a ] do
            M[j][h] := RowLengthen( galoisField,M[j][h],E[h][j],E[h][b] );            
        od;
    od;

    ## Step3 ##
    for k in [ 1 .. b ] do
        R[k][k] := D[k].vectors;
    od;
    for k_ in [ 1 .. b ] do
        k := b-k_+1;
        for j in [ 1 .. (k - 1) ] do
            tmp := CEX( galoisField,D[k].bitstring,B[j][k] );
            X := tmp[1];
            R[j][k] := tmp[2];
            if IsEmpty(X) then continue; fi;
            for l in [ k .. b ] do
                if not IsEmpty(R[k][l]) then
                    if IsEmpty(R[j][l]) then
                        R[j][l] :=X*R[k][l];
                    else
                        R[j][l] := R[j][l] + X*R[k][l];
                    fi;
                fi;
            od;
            for h in [ 1 .. a ] do
                if not IsEmpty(M[k][h]) then
                    if IsEmpty(M[j][h]) then
                        M[j][h] := X*M[k][h];
                    else
                        M[j][h] := M[j][h] + X*M[k][h];
                    fi;
                fi;
            od;
        od;
    od;
    Print("3\n");
    ###############################
    ###############################

    ## Write output
    # Begin with row-select bitstring named v
    v := [];
    rank := 0;
    w := 0*[1..(DimensionsMat(mat)[1]/a) ];
    for i in [ 1 .. a ] do
        if IsEmpty(E[i][b].rho) then
            tmp := w; 
        else
            tmp := E[i][b].rho;
        fi;
        v := Concatenation( v,tmp );
        rank := rank + Sum( tmp );
    od;

    #### Glueing the blocks of M
    B := NullMat( rank,Length(v),galoisField );
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
        w := DimensionsMat(mat)[1]/a;
        for i in [ 1 .. a ] do
            if IsEmpty(M[j][i]) then
                if IsEmpty(E[i][b].rho) then
                    tmp := w;
                else
                    tmp := Length( E[i][b].rho );
                fi;
                tmpC := tmpC + tmp; continue;
                #M[j][i] := NullMat( rows[j],tmp,galoisField );
            else
                M[j][i] := TransposedMat( RRF( galoisField,
                NullMat( Length(E[i][b].rho)-
                DimensionsMat(M[j][i])[2],
                DimensionsMat(M[j][i])[1],galoisField ),
                TransposedMat(M[j][i]),E[i][b].rho ) );
            fi;

            B{[tmpR .. tmpR + DimensionsMat(M[j][i])[1]-1 ]}{[tmpC .. tmpC + DimensionsMat(M[j][i])[2]-1 ]}
            := M[j][i];
            tmpC := tmpC + DimensionsMat(M[j][i])[2];
        od;
        tmpR := tmpR + rows[j];
    od;
   
    ### Glueing the blocks of R
    C := NullMat( rank,ncols-rank,galoisField );
    rows := [];
    w := [];
    for i in [ 1 .. b ] do
         rows[i] := 0;
         if IsEmpty(D[i].bitstring) then
             w := Concatenation( w,0*[1..nrows[i]] );
         else
             w := Concatenation( w,D[i].bitstring );
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
                 if not IsEmpty(D[j].bitstring) then
                     tmpC := tmpC + Sum( 1 - D[j].bitstring );
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

     C := TransposedMat( RRF( galoisField,TransposedMat(C), -IdentityMat( rank,galoisField ),w  ) );

    ## Glueing the blocks of K
    D := NullMat( Length(v)-rank,Length(v),galoisField );
    X := IdentityMat( Length(v)-rank,galoisField );
    rows := [];
    for j in [ 1 .. a ] do
        rows[j] := 0;
        for i in [ 1 .. a ] do
            if IsEmpty(K[j][i]) then
                rows[j] := Length(E[j][b].rho) - Sum(E[j][b].rho);
            else
                rows[j] := DimensionsMat(K[j][i])[1];
                break;
            fi;
        od;
    od;
    tmpR := 1;
    for j in [ 1 .. a ] do
        if rows[j]=0 then
            continue;
        fi;
        tmpC := 1;
        for i in [ 1 .. a ] do
            if IsEmpty(K[j][i]) then
                if IsEmpty(E[i][b].rho) then
                    #tmp := b;
                    tmp := DimensionsMat(mat)[1]/a;
                else
                    tmp := Length( E[i][b].rho );
                fi;
                tmpC := tmpC + tmp; continue;
                #K[j][i] := NullMat( rows[j],tmp,galoisField );
            else
                K[j][i] := TransposedMat( RRF( galoisField,
                    NullMat( Length(E[i][b].rho)-
                    DimensionsMat(K[j][i])[2],
                    DimensionsMat(K[j][i])[1],galoisField ),
                    TransposedMat(K[j][i]),E[i][b].rho ) );
            fi;

            D{[tmpR .. tmpR + DimensionsMat(K[j][i])[1]-1 ]}
            {[tmpC .. tmpC + DimensionsMat(K[j][i])[2]-1 ]}
            := K[j][i];
            tmpC := tmpC + DimensionsMat(K[j][i])[2];
        od;
        tmpR := tmpR + rows[j];
    od;

     D := TransposedMat( RRF( galoisField,X,
        TransposedMat( CEX( galoisField,v,D )[1] ),v ) );

    return rec( coeffs:=B,vectors:=C,relations:=D,
                pivotrows:=v,pivotcols:=w,rank:=rank);

end;

