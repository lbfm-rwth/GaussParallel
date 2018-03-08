# Contains a version of the elimination alg. for HPCGAP computing RREF and 
# a transformation, where the second step of the
#                    algorithm runs in parallel ( using HPCGAP's task arch. )


ClearUp := function( R,X,R_ )
    if IsEmpty(R_) or IsEmpty(X) then return R; fi;
    if IsEmpty(R) then
        return X*R_;
    else
        return R + X*R_;
    fi;
end;

ChiefParallelClearUp := function( galoisField,mat,a,b )
    ## inputs: a finite field, a matrix, natural nr. a,b to treat mat as axb matirx
    local   TaskListPreClearUp,
            TaskListClearUpM,
            TaskListClearUpR,
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
    dummyTask := RunTask( function() return []; end );
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

    ncols := DimensionsMat( mat )[2];
    C := ChopMatrix( galoisField,mat,a,b );
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
        D[k] := rec( remnant:=[],bitstring:=[] );
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
            tmp := ClearDown( galoisField,C[i][j],D[j],i );
            D[j] := tmp.D;
            A[i][j] := tmp.A;
            if j=1 then
                bs := rec( rho:=[],delta:=[],nr:=0 );
            else
                bs := E[i][j-1];
            fi;
            E[i][j] := Extend( A[i][j],bs,j );

            for k in [ j+1 .. b ] do
                tmp := UpdateRow(  galoisField,A[i][j],C[i][k],
                                    B[j][k],i );
                C[i][k] := tmp.C;
                B[j][k] := tmp.B;
            od;
            for h in [ 1 .. i ] do
                tmp := UpdateRowTrafo(  galoisField,A[i][j],K[i][h],
                                    M[j][h],E[h][j],i,h,j );
                K[i][h] := tmp.K;
                M[j][h] := tmp.M;
            od;
        od;    
    od;

    ## Step 2 ##
    for j in [ 1 .. b ] do
        for h in [ 1 .. a ] do
            M[j][h] := RowLengthen( galoisField,M[j][h],E[h][j],E[h][b] );            
        od;
    od;

    ## Step3 // parallel version ##
    for k in [ 1 .. b ] do
        R[k][k] := D[k].remnant;
    od;
    for k_ in [ 1 .. b ] do
        k := b-k_+1;
        for j in [ 1 .. (k - 1) ] do
            TaskListPreClearUp[j][k] := RunTask(
                CEX,
                galoisField,
                D[k].bitstring,
                B[j][k]
            );
            for l in [ k .. b ] do
                if l-k = 0 then
                    TaskListClearUpR[j][l][1] := ScheduleTask(
                        [   TaskListPreClearUp[j][l],
                            TaskListPreClearUp[j][k]
                        ],
                        ClearUp,
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
                        ClearUp,
                        TaskResult( TaskListClearUpR[j][l][l-k] ),
                        TaskResult( TaskListPreClearUp[j][k] )[1],
                        TaskResult( TaskListClearUpR[k][l][l-k] )
                    );
                fi;
            od;
            for h in [ 1 .. a ] do
                if k_ = 1 then
                    TaskListClearUpM[j][h][1] := ScheduleTask(
                        [
                            TaskListPreClearUp[j][k]
                        ],
                        ClearUp,
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
                        ClearUp,
                        TaskResult( TaskListClearUpM[j][h][k_-1] ),
                        TaskResult( TaskListPreClearUp[j][k] )[1],
                        TaskResult( TaskListClearUpM[k][h][k_-1] )
                    );
                fi;
            od;
        od;
    od;
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
            B{[tmpR .. tmpR + DimensionsMat(M[j][i])[1]-1 ]}
            {[tmpC .. tmpC + DimensionsMat(M[j][i])[2]-1 ]}
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

            C{[tmpR .. tmpR + DimensionsMat(R[i][j])[1]-1 ]}
             {[tmpC .. tmpC + DimensionsMat(R[i][j])[2]-1 ]}
             := R[i][j];
             tmpC := tmpC + DimensionsMat(R[i][j])[2];
        od;
        tmpR := tmpR + rows[i];
    od;    
    C := TransposedMat( 
        RRF( galoisField,TransposedMat(C), 
        -IdentityMat( rank,galoisField ),w  ) 
    );

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

    return rec( transformation:=B,remnant:=C,relations:=D,
                pivotrows:=v,pivotcols:=w,rank:=rank);

end;
