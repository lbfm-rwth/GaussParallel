# Install global functions of this package
# EchelonMat{, Transformation}Blockwise are wrappers for the function that does
# the actual work:
# DoEchelonMatTransformationBlockwise
InstallGlobalFunction(EchelonMatTransformationBlockwise,
function (mat)
    local result;
    result := DoEchelonMatTransformationBlockwise(mat, rec());
    return rec(
        vectors := result.vectors,
        heads := result.heads,
        coeffs := result.coeffs,
        relations := result.relations
    );
end
);

InstallGlobalFunction(EchelonMatBlockwise,
function (mat)
    local result;
    result := DoEchelonMatTransformationBlockwise(mat, rec());
    return rec(
        vectors := result.vectors,
        heads := result.heads
    );
end
);

InstallGlobalFunction(DoEchelonMatTransformationBlockwise,
function (mat, options)
    # options is a record that can specify
    # galoisField, IsHPC(removed for now), numberBlocksHeight, numberBlocksWidth,
    # withTrafo and verify
    local recnames, dim, numberBlocksHeight, numberBlocksWidth, a, b,
    withTrafo, galoisField, verify, C, ncols,
    TaskListClearDown, TaskListE, TaskListUpdateR, TaskListUpdateM,
    TaskListPreClearUp, TaskListClearUpR, TaskListClearUpM,
    A, B, D, E, K, M, R, X,
    nrows, ClearDownDeps, ExtendDeps, UpdateRowDeps, UpdateRowTrafoDeps,
    k, result, i, h, j, k_, l;

    recnames := Set(RecNames(options));

    dim := DimensionsMat(mat);
    Info(InfoGauss, 4, "Input checks");
    if ("numberBlocksHeight" in recnames)
            and ("numberBlocksWidth" in recnames) then
        numberBlocksHeight := options.numberBlocksHeight;
        numberBlocksWidth := options.numberBlocksWidth;
    else
        numberBlocksHeight := GAUSS_calculateBlocks(dim[1]);
        numberBlocksWidth := GAUSS_calculateBlocks(dim[2]);
    fi;
    #Hack
    a := numberBlocksHeight;
    b := numberBlocksWidth;

    if "withTrafo" in recnames then
        withTrafo := options.withTrafo;
    else
        withTrafo := true;
    fi;

    if "galoisField" in recnames then
        galoisField := options.galoisField;
    else
        galoisField := DefaultFieldOfMatrix(mat);
        if galoisField = fail then
            Error("Please specify the field of the matrix using the options parameter.");
        fi;
    fi;

    # verify only makes sense if withTrafo is true
    if "verify" in recnames then
        verify := options.verify;
    else
        verify := withTrafo;
    fi;
    if verify and not withTrafo then
        Error("can't verify the calculation without computing the ",
              "transformation matrix");
    fi;

    if not (HasIsField(galoisField) and IsField(galoisField)) then
        ErrorNoReturn("Wrong argument: The first parameter is not a field.");
    fi;
    if not IsMatrix(mat) then
        ErrorNoReturn("Wrong argument: The second parameter is not a matrix.");
    fi;
    if not (a in NonnegativeIntegers and b in NonnegativeIntegers) then
        ErrorNoReturn("Wrong argument: The third or fourth parameter is not a",
                      " nonnegative integer.");
    fi;

    ##Not supported yet
    if (DimensionsMat(mat)[1] mod a <> 0) then
        ErrorNoReturn("Variable: 'a' must divide number of rows");
    fi;
    if (DimensionsMat(mat)[2] mod b <> 0) then
        ErrorNoReturn("Variable: 'b' must divide number of columns");
    fi;

    Info(InfoGauss, 1, "The matrix is split into ", numberBlocksHeight,
        " blocks vertically and ", numberBlocksWidth, " horizontally.");

    Info(InfoGauss, 2, "------------ Start Chief ------------");
    Info(InfoGauss, 2, "Preparation");

    ##Preparation: Init and chopping the matrix mat

    C := GAUSS_ChopMatrix(galoisField, mat, a, b);
    ncols := DimensionsMat(mat)[2];

    TaskListClearDown := List([1 .. a], x -> List([1 .. b]));
    TaskListE := List([1 .. a], x -> List([1 .. b]));
    TaskListUpdateR := List(
        [1 .. a],
        x -> List(
            [1 .. b],
            x -> List(
                [1 .. b],
                x -> RunTask(function() return rec(C := [], B:=[]); end)
            )
        )
    );
    TaskListUpdateM := List(
        [1 .. a],
        x -> List(
            [1 .. b],
            x -> List(
                [1 .. a],
                x -> RunTask(function() return rec(M := [], K:=[]); end)
            )
        )
    );
    TaskListPreClearUp := List(
        [1..b],
        x -> []
    );
    TaskListClearUpR := List(
        [1..b],
        x -> List(
            [1..b],
            x -> []
        )
    );
    TaskListClearUpM := List(
        [1..b],
        x -> List(
            [1..a],
            x -> []
        )
    );

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
    for i in [1 .. a] do
        A[i] := FixedAtomicList(b, 0);
        E[i] := FixedAtomicList(b, 0);
        K[i] := FixedAtomicList(a, 0);
        for k in [1 .. b] do
            A[i][k] := MakeReadOnlyOrImmutableObj(
                rec(A := [], M := [], E := [], K := [],
                rho := [], lambda := [])
            );
            E[i][k] := MakeReadOnlyOrImmutableObj(
                rec(rho:=[], delta:=[], nr:=0)
            );
        od;
        for h in [1 .. a] do
            K[i][h] := MakeReadOnlyOrImmutableObj([]);
        od;
    od;
    for i in [1 .. b] do
        M[i] := FixedAtomicList(a, 0);
        for k in [1 .. a] do
            M[i][k] := MakeReadOnlyOrImmutableObj([]);
        od;
    od;
    for k in [1 .. b] do
        D[k] := MakeReadOnlyOrImmutableObj(
            rec(vectors := [], bitstring := [])
        );
        B[k] := FixedAtomicList(b, 0);
        R[k] := FixedAtomicList(b, 0);
        for j in [1 .. b] do
            B[k][j] := MakeReadOnlyOrImmutableObj([]);
            R[k][j] := MakeReadOnlyOrImmutableObj([]);
        od;
    od;
    # X is only used by the main thread so we don't need to make it threadsafe
    # b x b
    X := [];
    for i in [1 .. b] do
        X[k] := [];
        for j in [1 .. b] do
            X[k][j] := [];
        od;
    od;
    # nrows is only used when glueing together R
    nrows := [];
    for i in [1 .. b] do
        nrows[i] := DimensionsMat(C[1][i])[2];
    od;
    ###############################
    ###############################

    ## Step 1 ##
    Info(InfoGauss, 2, "Step 1");
    for i in [1 .. a] do
        for j in [1 .. b] do
            Info(InfoGauss, 3, "ClearDownParameters ", i, " ", j);
            ClearDownDeps := GAUSS_ClearDownDependencies(
                i, j, TaskListClearDown, TaskListUpdateR
            );
            TaskListClearDown[i][j] := ScheduleTask(
                ClearDownDeps, GAUSS_ClearDown_destructive,
                galoisField, C, D, A, i, j
            );
            Info(InfoGauss, 3, "ExtendParameters ", i, " ", j);
            ExtendDeps := GAUSS_ExtendDependencies(
                i, j, TaskListClearDown, TaskListE
            );
            TaskListE[i][j] := ScheduleTask(
                ExtendDeps, GAUSS_Extend_destructive, A, E, i, j
            );
            Info(InfoGauss, 3, "UpdateRowParameters ", i, " ", j);
            for k in [j+1 .. b] do
                UpdateRowDeps := GAUSS_UpdateRowDependencies(
                    i, j, k, TaskListClearDown, TaskListUpdateR
                );
                TaskListUpdateR[i][j][k] := ScheduleTask(
                    UpdateRowDeps, GAUSS_UpdateRow_destructive,
                    galoisField, A, C, B, i, j, k
                );
            od;

            if withTrafo then
                Info(InfoGauss, 3, "UpdateRowTrafoParameters ", i, " ", j);
                for h in [1 .. i] do
                    UpdateRowTrafoDeps := GAUSS_UpdateRowTrafoDependencies(
                        i, j, h, TaskListClearDown, TaskListE, TaskListUpdateM
                    );
                    TaskListUpdateM[i][j][h] := ScheduleTask(
                        UpdateRowTrafoDeps, GAUSS_UpdateRowTrafo_destructive,
                        galoisField, A, K, M, E, i, h, j
                    );
                od;
            fi;
        od;
    od;

    Info(InfoGauss, 3, "Before WaitTask");
    WaitTask(Concatenation(TaskListClearDown));
    WaitTask(Concatenation(TaskListE));
    WaitTask(Concatenation(List(TaskListUpdateR, Concatenation)));
    WaitTask(Concatenation(List(TaskListUpdateM, Concatenation)));
    if  withTrafo then
        ## Step 2 ##
        Info(InfoGauss, 2, "Step 2");
        for j in [1 .. b] do
            for h in [1 .. a] do
                M[j][h] := MakeReadOnlyOrImmutableObj(
                    GAUSS_RowLengthen(galoisField, M[j][h], E[h][j], E[h][b])
                );
            od;
        od;
    fi;

    ## Step3 ##
    Info(InfoGauss, 2, "Step 3");
    for k in [1 .. b] do
        R[k][k] := ShallowCopy(D[k].vectors);
        MakeReadOnlyOrImmutableObj(R[k][k]);
    od;
    Info(InfoGauss, 2, "ClearUpR and possibly ClearUpM");
    for k_ in [1 .. b] do
        k := b - k_ + 1;
        for j in [1 .. (k - 1)] do
            TaskListPreClearUp[j][k] := RunTask(
                GAUSS_PreClearUp,
                R, galoisField, D, B, j, k
            );
            for l in [k .. b] do
                Info(InfoGauss, 3,
                     "ClearUpR j, k, l: ", j, " ", k, " ", l);
                if l-k = 0 then
                    TaskListClearUpR[j][l][1] := ScheduleTask(
                        [TaskListPreClearUp[j][k]],
                        GAUSS_ClearUp_destructive,
                        R, TaskResult(TaskListPreClearUp[j][k]), j, k, l
                    );
                else
                    TaskListClearUpR[j][l][l-k+1] := ScheduleTask(
                        [
                            TaskListClearUpR[j][l][l-k],
                            TaskListClearUpR[k][l][l-k],
                            TaskListPreClearUp[j][k]
                        ],
                        GAUSS_ClearUp_destructive,
                        R, TaskResult(TaskListPreClearUp[j][k]), j, k, l
                    );
                fi;
            od;
            if  withTrafo then
                for h in [1 .. a] do
                    Info(InfoGauss, 3,
                         "ClearUpM j, k, h: ", j, " ", k, " ", h);
                    if k_ = 1 then
                        TaskListClearUpM[j][h][1] := ScheduleTask(
                            [TaskListPreClearUp[j][k]],
                            GAUSS_ClearUp_destructive,
                            M, TaskResult(TaskListPreClearUp[j][k]), j, k, h
                        );
                    else
                        TaskListClearUpM[j][h][k_] := ScheduleTask(
                            [
                                TaskListClearUpM[j][h][k_-1],
                                TaskListClearUpM[k][h][k_-1],
                                TaskListPreClearUp[j][k]
                            ],
                            GAUSS_ClearUp_destructive,
                            M, TaskResult(TaskListPreClearUp[j][k]), j, k, h
                        );
                    fi;
                od;
            fi;
        od;
    od;

    # Wait for the computations of step 3.
    WaitTask(Concatenation(TaskListPreClearUp));
    WaitTask(Concatenation(List(TaskListClearUpR, Concatenation)));
    if  withTrafo then
        WaitTask(Concatenation(List(TaskListClearUpM, Concatenation)));
    fi;
    # Computations are finished. Now prepare the output.
    result := GAUSS_WriteOutput(mat, a, b, ncols, nrows, galoisField, D, R, M,
                                E, K, withTrafo);
    if verify and not IsMatrixInRREF(result.transformation * mat) then
        Error("Result verification failed! Result is not in RREF!\n",
              "Please report this bug to\n",
              "https://github.com/lbfm-rwth/GaussPar#contact");
    fi;
    return result;
end);
