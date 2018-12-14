# Contains a version of the elimination algorithm for HPCGAP computing RREF and a 
# transformation, running completely in parallel.


GAUSS_ClearUp := function( R,X,R_ )
    if IsEmpty(R_) or IsEmpty(X) then return R; fi;
    if IsEmpty(R) then
        return X*R_;
    else
        return R + X*R_;
    fi;
end;

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

Chief := function( galoisField,mat,a,b,IsHPC )
    ## inputs: a finite field, a matrix, natural nr. a,b to treat mat as axb matirx
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
        TaskListPreClearUp := List(
            [ 1..b ],
            x -> List( [ 1..b ], x -> dummyTask )
        );
        TaskListClearUpR := List(
            [ 1..b ],
            x -> List( 
                [ 1..b ], 
                x -> List( [ 1..b ], x -> dummyTask )
            )
        );
        TaskListClearUpM := List(
            [ 1..b ],
            x -> List( 
                [ 1..a ], 
                x -> List( [ 1..b ], x -> dummyTask )
            )
        );
	fi;
        
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
    Info(InfoGauss, 2, "Step 1");
    for i in [ 1 .. a ] do
        for j in [ 1 .. b ] do
		    if IsHPC then
          		Info(InfoGauss, 3, "ClearDownParameters ", i, " ", j);
			    ClearDownInput := GAUSS_ClearDownParameters(i, j, C, TaskListClearDown,
				    TaskListUpdateR, galoisField);
			    TaskListClearDown[i][j] := ScheduleTask(
				    ClearDownInput[1],
				    GAUSS_ClearDown,
				    ClearDownInput[2],
				    ClearDownInput[3],
				    ClearDownInput[4],
				    ClearDownInput[5]
			    );
    
	        	    Info(InfoGauss, 3, "ExtendParameters ", i, " ", j);
			    ExtendInput := GAUSS_ExtendParameters(i, j, TaskListClearDown, TaskListE);
			    TaskListE[i][j] := ScheduleTask(
				    ExtendInput[1],
				    GAUSS_Extend,
				    ExtendInput[2],
				    ExtendInput[3],
				    ExtendInput[4]
			    );

				Info(InfoGauss, 3, "UpdateRowParameters ", i, " ", j);
		    else
		        tmp := GAUSS_ClearDown( galoisField,C[i][j],D[j],i );
		        D[j] := tmp.D;
		        A[i][j] := tmp.A;
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
				        UpdateRowInput[1],
				        GAUSS_UpdateRow,
				        UpdateRowInput[2],
				        UpdateRowInput[3],
				        UpdateRowInput[4],
				        UpdateRowInput[5],
				        UpdateRowInput[6]
			        );
		        else
		                tmp := GAUSS_UpdateRow(  galoisField,A[i][j],C[i][k],
		                                    B[j][k],i );
		                C[i][k] := tmp.C;
		                B[j][k] := tmp.B;
		        fi;
            od;

            Info(InfoGauss, 3, "UpdateRowTrafoParameters ", i, " ", j);
            for h in [ 1 .. i ] do
		        if IsHPC then
			        UpdateRowTrafoInput := GAUSS_UpdateRowTrafoParameters(i, j, h, TaskListClearDown, TaskListE, TaskListUpdateM, galoisField);
            	    TaskListUpdateM[i][j][h] := ScheduleTask(
					    UpdateRowTrafoInput[1],
                	    GAUSS_UpdateRowTrafo,
					    UpdateRowTrafoInput[2],
					    UpdateRowTrafoInput[3],
					    UpdateRowTrafoInput[4],
					    UpdateRowTrafoInput[5],
					    UpdateRowTrafoInput[6],
					    UpdateRowTrafoInput[7],
					    UpdateRowTrafoInput[8],
					    UpdateRowTrafoInput[9]
				    ); 
		        else
		            tmp := GAUSS_UpdateRowTrafo(  galoisField,A[i][j],K[i][h],
		                                M[j][h],E[h][j],i,h,j );
		            K[i][h] := tmp.K;
		            M[j][h] := tmp.M;
		        fi;
            od;
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
                for h in [ 1 .. i ] do
                    K[i][h] := TaskResult( TaskListUpdateM[i][j][h] ).K;
                    M[j][h] := TaskResult( TaskListUpdateM[i][j][h] ).M;
                od;
            od;
        od;
	fi;

    ## Step 2 ##
    Info(InfoGauss, 2, "Step 2");
    for j in [ 1 .. b ] do
        for h in [ 1 .. a ] do
            M[j][h] := GAUSS_RowLengthen( galoisField,M[j][h],E[h][j],E[h][b] );            
        od;
    od;

    ## Step3 ##
    for k in [ 1 .. b ] do
        R[k][k] := D[k].vectors;
    od;
    for k_ in [ 1 .. b ] do
        k := b-k_+1;
        for j in [ 1 .. (k - 1) ] do
		    if IsHPC then
                TaskListPreClearUp[j][k] := RunTask(
                    GAUSS_CEX,
                    galoisField,
                    D[k].bitstring,
                    B[j][k]
                );
		    else
                tmp := GAUSS_CEX( galoisField,D[k].bitstring,B[j][k] );
                X := tmp[1];
                R[j][k] := tmp[2];
                if IsEmpty(X) then continue; fi;
		    fi;

            for l in [ k .. b ] do
		        if IsHPC then
                        if l-k = 0 then
                            TaskListClearUpR[j][l][1] := ScheduleTask(
                                [   TaskListPreClearUp[j][l],
                                    TaskListPreClearUp[j][k]
                                ],
                                GAUSS_ClearUp,
                                TaskResult( TaskListPreClearUp[j][l] )[2],
                                TaskResult( TaskListPreClearUp[j][k] )[1],
                                R[k][l]
                            ); 
                        else
                            TaskListClearUpR[j][l][l-k+1] := ScheduleTask(
                                [   TaskListClearUpR[j][l][l-k],
                                    TaskListClearUpR[k][l][l-k],
                                    TaskListPreClearUp[j][k]
                                ],
                                GAUSS_ClearUp,
                                TaskResult( TaskListClearUpR[j][l][l-k] ),
                                TaskResult( TaskListPreClearUp[j][k] )[1],
                                TaskResult( TaskListClearUpR[k][l][l-k] )
                            );
                        fi;
		        else
                    if not IsEmpty(R[k][l]) then
                        if IsEmpty(R[j][l]) then
                            R[j][l] :=X*R[k][l];
                        else
                            R[j][l] := R[j][l] + X*R[k][l];
                        fi;
                    fi;
		        fi;
            od;

            for h in [ 1 .. a ] do
		        if IsHPC then
                        if k_ = 1 then
                            TaskListClearUpM[j][h][1] := ScheduleTask(
                                [
                                    TaskListPreClearUp[j][k]
                                ],
                                GAUSS_ClearUp,
                                M[j][h],
                                TaskResult( TaskListPreClearUp[j][k] )[1],
                                M[k][h]
                            ); 
                        else
                            TaskListClearUpM[j][h][k_] := ScheduleTask(
                                [   TaskListClearUpM[j][h][k_-1],
                                    TaskListClearUpM[k][h][k_-1],
                                    TaskListPreClearUp[j][k]
                                ],
                                GAUSS_ClearUp,
                                TaskResult( TaskListClearUpM[j][h][k_-1] ),
                                TaskResult( TaskListPreClearUp[j][k] )[1],
                                TaskResult( TaskListClearUpM[k][h][k_-1] )
                            );
                        fi;
		        else
                    if not IsEmpty(M[k][h]) then
                        if IsEmpty(M[j][h]) then
                            M[j][h] := X*M[k][h];
                        else
                            M[j][h] := M[j][h] + X*M[k][h];
                        fi;
                    fi;
		        fi;
            od;
        od;
    od;
	
	if IsHPC then
        WaitTask( Concatenation( TaskListPreClearUp ) );
        WaitTask( Concatenation( List( TaskListClearUpR,Concatenation ) ) );
        WaitTask( Concatenation( List( TaskListClearUpM,Concatenation ) ) );
        for i in [ 1 .. b ] do
            k := b - i + 1;
            for j in [ 1 .. k-1 ] do
                for h in [ k .. b ] do
                    R[j][h] := TaskResult( TaskListClearUpR[j][h][h-k+1] );
                od;
                for h in [ 1 .. a ] do
                    M[j][h] := TaskResult( TaskListClearUpM[j][h][i] );
                od;
            od;
        od;
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

    heads := GAUSS_createHeads(v, w, DimensionsMat(mat)[2]);

    return rec( coeffs:=B,vectors:=C,relations:=D,
                pivotrows:=v,pivotcols:=w,rank:=rank,
                heads:=heads);

end;

