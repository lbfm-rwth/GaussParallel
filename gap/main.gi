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
    local recnames, recognisedOptions, nrows, ncols, numberBlocksHeight,
    numberBlocksWidth, a, b, withTrafo, galoisField, verify, tmp, C,
    nrRowsPerBlockRow, nrColsPerBlockCol, copyMat,
    TaskListClearDown, TaskListE, TaskListUpdateR, TaskListUpdateM,
    TaskListPreClearUp, TaskListClearUpR, TaskListClearUpM,
    A, B, D, E, K, M, R, X,
    ClearDownDeps, ExtendDeps, UpdateRowDeps, UpdateRowTrafoDeps,
    k, result, i, h, j, k_, l, ClearUpDeps, isChopped;

    recnames := Set(RecNames(options));

    Info(InfoGauss, 4, "Input checks");
    if ("numberBlocks" in recnames) then
	if ("numberBlocksHeight" in recnames) or ("numberBlocksWidth" in recnames) then
	    ErrorNoReturn("Conflicting parameters: Can't use numberBlocks and numberBlocksWidth/numberBlocksHeight simultaneously");
	fi;
        numberBlocksHeight := options.numberBlocks;
        numberBlocksWidth := options.numberBlocks;		
    elif ("numberBlocksHeight" in recnames)
            and ("numberBlocksWidth" in recnames) then
        numberBlocksHeight := options.numberBlocksHeight;
        numberBlocksWidth := options.numberBlocksWidth;
    else
        numberBlocksHeight := 1; 
        numberBlocksWidth := 1;
    fi;
    #Hack
    a := numberBlocksHeight;
    b := numberBlocksWidth;

    if "withTrafo" in recnames then
        withTrafo := options.withTrafo;
    else
        withTrafo := true;
    fi;
    
    if "isChopped" in recnames then
        isChopped := options.isChopped;
    else
        isChopped := false;
    fi;
    if not isChopped and not IsMatrix(mat) then
        ErrorNoReturn("<mat> is not a matrix.");
    fi;
    if isChopped then
	if not (("numberBlocks" in recnames) or ("numberBlocksHeight" in recnames)) then
			ErrorNoReturn("Insufficient Paramters: Can't use isChopped option without specifying numberBlocks or both of numberBlocksHeight and numberBlocksWidth."); 
	fi;
        nrows := 0;
        ncols := 0;
        for i in [ 1 .. numberBlocksHeight ] do
            nrows := nrows + NrRows(mat[i][1]);
        od;     
        for j in [ 1 .. numberBlocksWidth ] do
            ncols := ncols + NrCols(mat[1][j]);
        od;
    else
        nrows := NrRows(mat);
        ncols := NrCols(mat);
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
        verify := false;
    fi;
    if verify and isChopped then
        Error("Can't use verify and isChopped simultaneously");
    fi;
    if verify and not withTrafo then
        Error("can't verify the calculation without computing the ",
              "transformation matrix");
    fi;

    if not (HasIsField(galoisField) and IsField(galoisField)) then
        ErrorNoReturn("<galoisField> is not a field.");
    fi;
    if not (a in NonnegativeIntegers and b in NonnegativeIntegers) then
        ErrorNoReturn("<numberBlocksWidth> and <numberBlocksHeight> must be ",
                      " nonnegative integers.");
    fi;
    if not (numberBlocksHeight <= nrows and numberBlocksWidth <= ncols) then
        ErrorNoReturn("<numberBlocksHeight> and <numberBlocksWidth> must be",
                      " less or equal than the number of rows and columns",
                      " respectively");
    fi;
    # Check for invalid components of `options`. Note that IsHPC is a dummy
    # value that is still used by the tests but doesn't do anything.
    recognisedOptions := ["numberBlocksHeight", "numberBlocksWidth",
                          "withTrafo", "galoisField", "verify",
                          "IsHPC"];
    if not IsSubset(recognisedOptions, recnames) then
        Print("Warning: unrecognised options: ",
              Difference(recnames, recognisedOptions), "\n");
    fi;

    Info(InfoGauss, 2, "The matrix is split into ", numberBlocksHeight,
        " blocks vertically and ", numberBlocksWidth, " horizontally.");

    Info(InfoGauss, 2, "------------ Start Chief ------------");
    Info(InfoGauss, 2, "Preparation");

    ##Preparation: Init and chopping the matrix mat
    if isChopped then
        copyMat := FixedAtomicList(a,0);
        for i in [ 1 .. a ] do
            copyMat[i] := FixedAtomicList(b,0);
            for j in [ 1 .. b ] do
                copyMat[i][j] := MutableCopyMat(mat[i][j]);
                ConvertToMatrixRepNC(copyMat[i][j], galoisField);
            od;
        od;
    else 
        copyMat := MutableCopyMat(mat);
        ConvertToMatrixRepNC(copyMat, galoisField); 
    fi;
    tmp := GAUSS_ChopMatrix(galoisField, copyMat, a, b, isChopped);
    C := tmp.mat;
    nrRowsPerBlockRow := tmp.rowsList;
    nrColsPerBlockCol := tmp.colsList;

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
        nrows[i] := NrCols(C[1][i]);
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
                ClearUpDeps := GAUSS_ClearUpDependencies(
                    j, k, l, l-k, TaskListPreClearUp, TaskListClearUpR
                );                               
                TaskListClearUpR[j][l][l-k+1] := ScheduleTask(
                    ClearUpDeps, GAUSS_ClearUp_destructive,
                    R, TaskResult(TaskListPreClearUp[j][k]), j, k, l
                );
            od;
            if  withTrafo then
                for h in [1 .. a] do
                    Info(InfoGauss, 3,
                         "ClearUpM j, k, h: ", j, " ", k, " ", h);
                    ClearUpDeps := GAUSS_ClearUpDependencies(
                        j, k, h, k_-1, TaskListPreClearUp, TaskListClearUpM
                    );                               
                    TaskListClearUpM[j][h][k_] := ScheduleTask(
                        ClearUpDeps, GAUSS_ClearUp_destructive,
                        M, TaskResult(TaskListPreClearUp[j][k]), j, k, h
                    );
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
    result := GAUSS_WriteOutput(mat, a, b, nrColsPerBlockCol, nrRowsPerBlockRow, galoisField, D, R, M,
                                E, K, withTrafo);
	
    #the transformation matrix is just Concatenation(result.coeffs, result.relations);
    if verify and not IsMatrixInRREF(Concatenation(result.coeffs, result.relations) * mat) then
        Error("Result verification failed! Result is not in RREF!\n",
              "Please report this bug to\n",
              "https://github.com/lbfm-rwth/GaussPar#contact");
    fi;
    return result;
end);
