# Contains a version of the elimination algorithm for HPCGAP computing RREF and a
# transformation, running completely in parallel.


GAUSS_createHeads := function( pivotrows, pivotcols, width )
    # inputs: list that contains the row numbers of the pivot rows and
    # list that contains the column numbers of the pivot cols and
    # width of matrix
    local result, i, currentPivot;

    if not Length(pivotrows) = Length(pivotcols) then
        return [];
    fi;

    currentPivot := 0;
    result := ListWithIdenticalEntries( width, 0 );

    for i in [1 .. width] do
        if i in pivotcols then
            currentPivot := currentPivot + 1;
            result[i] := pivotrows[currentPivot];
        else
            result[i] := 0;
        fi;
    od;

    return result;
end;

Chief := function( galoisField,mat,a,b,IsHPC,withTrafo )
    ## inputs: a finite field, a matrix, number of vertical blocks, number of horizontal blocks
    local   TaskListPreClearUp,
            TaskListClearDown,
            TaskListClearUpR,
            TaskListClearUpM,
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
            h,
            ClearDownInput,
            ExtendInput,
            UpdateRowInput,
            UpdateRowTrafoInput,
            heads;

    ##Preparation: Init and chopping the matrix mat
    Info(InfoGauss, 2, "------------ Start Chief ------------");
    Info(InfoGauss, 2, "Preparation");

    Info(InfoGauss, 4, "Input checks");
    if not (HasIsField(galoisField) and IsField(galoisField)) then
        ErrorNoReturn("Wrong argument: The first parameter is not a field.");
    fi;
    if not IsMatrix(mat) then
        ErrorNoReturn("Wrong argument: The second parameter is not a matrix.");
    fi;
    if not (a in NonnegativeIntegers and b in NonnegativeIntegers) then
        ErrorNoReturn("Wrong argument: The third or fourth parameter is not a nonnegative integer.");
    fi;
    if not IsBool(IsHPC) then
        ErrorNoReturn("Wrong argument: The fifth parameter is not a boolean.");
    fi;

    ##Not supported yet
    if ( DimensionsMat(mat)[1] mod a <> 0) then
        ErrorNoReturn("Variable: 'a' must divide number of rows" );
    fi;
    if ( DimensionsMat(mat)[2] mod b <> 0) then
        ErrorNoReturn("Variable: 'b' must divide number of columns" );
    fi;

    C := GAUSS_ChopMatrix( galoisField,mat,a,b );
    ncols := DimensionsMat( mat )[2];

    if IsHPC then
        TaskListClearDown := List(
            [ 1 .. a ],
            x -> List( [ 1 .. b ] )
        );
        TaskListE := List(
            [ 1 .. a ],
            x -> List( [ 1 .. b ] )
        );
        TaskListUpdateR := List(
            [ 1 .. a ],
            x -> List(
                [ 1 .. b ],
                x -> List( [ 1 .. b ], x -> RunTask( function() return rec( C := [],B:=[] ); end ) )
            )
        );
        TaskListUpdateM := List(
            [ 1 .. a ],
            x -> List(
                [ 1 .. b ],
                x -> List( [ 1 .. a ], x -> RunTask( function() return rec( M := [],K:=[] ); end ) )
            )
        );
        TaskListPreClearUp := List(
            [ 1..b ],
            x -> []
        );
        TaskListClearUpR := List(
            [ 1..b ],
            x -> List(
                [ 1..b ],
                x -> []
            )
        );
        TaskListClearUpM := List(
            [ 1..b ],
            x -> List(
                [ 1..a ],
                x -> []
            )
        );
    fi;

    # List dimensions:
    # a x b
    A := FixedAtomicList(a, 0);
    # b x b
    B := FixedAtomicList(b, 0);
    # b
    D := FixedAtomicList(b, 0);
    # a x b
    E := FixedAtomicList(a, 0);
    # a x a
    K := FixedAtomicList(a, 0);
    # b x a
    M := FixedAtomicList(b, 0);
    # b x b
    R := FixedAtomicList(b, 0);
    for i in [ 1 .. a ] do
        A[i] := FixedAtomicList(b, 0);
        E[i] := FixedAtomicList(b, 0);
        K[i] := FixedAtomicList(a, 0);
        for k in [ 1 .. b ] do
            A[i][k] := MakeReadOnlyObj(MakeImmutable(
                rec(A := [], M := [], E := [], K := [],
                rho := [], lambda := [])
            ));
            E[i][k] := MakeReadOnlyObj(MakeImmutable(
                rec( rho:=[],delta:=[],nr:=0 )
            ));
        od;
        for h in [ 1 .. a ] do
            K[i][h] := MakeReadOnlyObj(MakeImmutable([]));
        od;
    od;
    for i in [ 1 .. b ] do
        M[i] := FixedAtomicList(a, 0);
        for k in [ 1 .. a ] do
            M[i][k] := MakeReadOnlyObj(MakeImmutable([]));
        od;
    od;
    for k in [ 1 .. b ] do
        D[k] := MakeReadOnlyObj(MakeImmutable(
            rec(vectors := [], bitstring := [])
        ));
        B[k] := FixedAtomicList(b, 0);
        R[k] := FixedAtomicList(b, 0);
        for j in [ 1 .. b ] do
            B[k][j] := MakeReadOnlyObj(MakeImmutable([]));
            R[k][j] := MakeReadOnlyObj(MakeImmutable([]));
        od;
    od;
    # X is only used by the main thread so we don't need to make it threadsafe
    # b x b
    X := [];
    for i in [ 1 .. b ] do
        X[k] := [];
        for j in [ 1 .. b ] do
            X[k][j] := [];
        od;
    od;
    # nrows is only used when glueing together R
    nrows := [];
    for i in [ 1 .. b ] do
        nrows[i] := DimensionsMat(C[1][i])[2];
    od;
    ###############################
    ###############################

    ## Step 1 ##
    Info(InfoGauss, 2, "Step 1");
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
            if IsHPC then
                  Info(InfoGauss, 3, "ClearDownParameters ", i, " ", j);
                ClearDownInput := GAUSS_ClearDownParameters(i, j, C, TaskListClearDown,
                    TaskListUpdateR, galoisField);
                TaskListClearDown[i][j] := ScheduleTask(
                    ClearDownInput.dependencies,
                    GAUSS_ClearDown,
                    ClearDownInput.parameters.galoisField,
                    ClearDownInput.parameters.C,
                    ClearDownInput.parameters.D,
                    ClearDownInput.parameters.i
                );

                    Info(InfoGauss, 3, "ExtendParameters ", i, " ", j);
                ExtendInput := GAUSS_ExtendParameters(i, j, TaskListClearDown, TaskListE);
                TaskListE[i][j] := ScheduleTask(
                    ExtendInput.dependencies,
                    GAUSS_Extend,
                    ExtendInput.parameters.A,
                    ExtendInput.parameters.E,
                    ExtendInput.parameters.flag
                );

                Info(InfoGauss, 3, "UpdateRowParameters ", i, " ", j);
            else
                GAUSS_ClearDown_destructive( galoisField,C,D,A,i,j );
                if j=1 then
                bs := rec( rho:=[],delta:=[],nr:=0 );
                else
                bs := E[i][j-1];
                fi;
                E[i][j] := GAUSS_Extend( A[i][j],bs,j );
            fi;

            for k in [ j+1 .. b ] do
                if IsHPC then
                    UpdateRowInput := GAUSS_UpdateRowParameters(i, j, k, C, TaskListClearDown,
                        TaskListUpdateR, galoisField);
                    TaskListUpdateR[i][j][k] := ScheduleTask(
                        UpdateRowInput.dependencies,
                        GAUSS_UpdateRow,
                        UpdateRowInput.parameters.galoisField,
                        UpdateRowInput.parameters.A,
                        UpdateRowInput.parameters.C,
                        UpdateRowInput.parameters.B,
                        UpdateRowInput.parameters.i
                    );
                else
                        GAUSS_UpdateRow_destructive(  galoisField,A,C,B,i,j,k );
                fi;
            od;

            if withTrafo then
                Info(InfoGauss, 3, "UpdateRowTrafoParameters ", i, " ", j);
                for h in [ 1 .. i ] do
                    if IsHPC then
                        UpdateRowTrafoInput := GAUSS_UpdateRowTrafoParameters(i, j, h, TaskListClearDown, TaskListE, TaskListUpdateM, galoisField);
                        TaskListUpdateM[i][j][h] := ScheduleTask(
                            UpdateRowTrafoInput.dependencies,
                            GAUSS_UpdateRowTrafo,
                            UpdateRowTrafoInput.parameters.galoisField,
                            UpdateRowTrafoInput.parameters.A,
                            UpdateRowTrafoInput.parameters.K,
                            UpdateRowTrafoInput.parameters.M,
                            UpdateRowTrafoInput.parameters.E,
                            UpdateRowTrafoInput.parameters.i,
                            UpdateRowTrafoInput.parameters.k,
                            UpdateRowTrafoInput.parameters.j
                        );
                    else
                        GAUSS_UpdateRowTrafo_destructive(  galoisField,A,K,M,E,i,h,j );
                    fi;
                od;
            fi;
        od;
    od;

    if IsHPC then
    Info(InfoGauss, 3, "Before WaitTask");
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
                if withTrafo then
                    for h in [ 1 .. i ] do
                        K[i][h] := TaskResult( TaskListUpdateM[i][j][h] ).K;
                        M[j][h] := TaskResult( TaskListUpdateM[i][j][h] ).M;
                    od;
                fi;
            od;
        od;
    fi;

    if  withTrafo then
        ## Step 2 ##
        Info(InfoGauss, 2, "Step 2");
        for j in [ 1 .. b ] do
            for h in [ 1 .. a ] do
                M[j][h] := GAUSS_RowLengthen( galoisField,M[j][h],E[h][j],E[h][b] );
            od;
        od;
    fi;

    ## Step3 ##
    Info(InfoGauss, 2, "Step 3");
    for k in [ 1 .. b ] do
        R[k][k] := ShallowCopy(D[k].vectors);
        MakeReadOnlyObj(R[k][k]); # FIXME: do we need to do this?
    od;
    Info(InfoGauss, 2, "CLearUpR");
    for k_ in [ 1 .. b ] do
        k := b-k_+1;
        for j in [ 1 .. (k - 1) ] do
            if IsHPC then
                TaskListPreClearUp[j][k] := RunTask(
                    GAUSS_PreClearUp,
                    R, galoisField, D, B, j, k
                );
            else
                X := GAUSS_PreClearUp( R,galoisField,D,B,j,k );
                if IsEmpty(X) then continue; fi;
            fi;

            for l in [ k .. b ] do
                if IsHPC then
                        if l-k = 0 then
                            TaskListClearUpR[j][l][1] := ScheduleTask(
                                [ TaskListPreClearUp[j][k] ],
                                GAUSS_ClearUp_destructive,
                                R,
                                TaskResult( TaskListPreClearUp[j][k] ),
                                j,
                                k,
                                l
                            );
                        else
                            TaskListClearUpR[j][l][l-k+1] := ScheduleTask(
                                [
                                    TaskListClearUpR[j][l][l-k],
                                    TaskListClearUpR[k][l][l-k],
                                    TaskListPreClearUp[j][k]
                                ],
                                GAUSS_ClearUp_destructive,
                                R,
                                TaskResult( TaskListPreClearUp[j][k] ),
                                j,
                                k,
                                l
                            );
                        fi;
                else
                    GAUSS_ClearUp_destructive( R,X,j,k,l );
                fi;
            od;
            
            if  withTrafo then
                Info(InfoGauss, 2, "CLearUpM");
                for h in [ 1 .. a ] do
                    if IsHPC then
                            if k_ = 1 then
                                TaskListClearUpM[j][h][1] := ScheduleTask(
                                    [
                                        TaskListPreClearUp[j][k]
                                    ],
                                    GAUSS_ClearUp_destructive,
                                    M,
                                    TaskResult( TaskListPreClearUp[j][k] ),
                                    j,
                                    k,
                                    h
                                );
                            else
                                TaskListClearUpM[j][h][k_] := ScheduleTask(
                                    [   TaskListClearUpM[j][h][k_-1],
                                        TaskListClearUpM[k][h][k_-1],
                                        TaskListPreClearUp[j][k]
                                    ],
                                    GAUSS_ClearUp_destructive,
                                    M,
                                    TaskResult( TaskListPreClearUp[j][k] ),
                                    j,
                                    k,
                                    h
                                );
                            fi;
                    else
                                GAUSS_ClearUp_destructive( M,X,j,k,h );
                    fi;
                od;
            fi;
        od;
    od;

    if IsHPC then
        WaitTask( Concatenation( TaskListPreClearUp ) );
        WaitTask( Concatenation( List( TaskListClearUpR,Concatenation ) ) );
        if  withTrafo then
            WaitTask( Concatenation( List( TaskListClearUpM,Concatenation ) ) );
        fi;
    fi;

    ###############################
    ###############################

    ## Write output
    Info(InfoGauss, 2, "Write output");
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

    if  withTrafo then
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
                    M[j][i] := TransposedMat( GAUSS_RRF( galoisField,
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
    fi;
    
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

     C := TransposedMat( GAUSS_RRF( galoisField,TransposedMat(C), -IdentityMat( rank,galoisField ),w  ) );

    if withTrafo then
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
                    K[j][i] := TransposedMat( GAUSS_RRF( galoisField,
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
        D := TransposedMat( GAUSS_RRF( galoisField,X,
            TransposedMat( GAUSS_CEX( galoisField,v,D )[1] ),v ) );
    fi;
    

    heads := GAUSS_createHeads(v, w, DimensionsMat(mat)[2]);
    if withTrafo then
        return rec( coeffs := -B, vectors := -C, relations := D,
                pivotrows := v, pivotcols := w, rank := rank,
                heads := heads );
        
    fi;
    return rec( vectors := -C,pivotrows := v, pivotcols := w, rank := rank, heads := heads );
end;

