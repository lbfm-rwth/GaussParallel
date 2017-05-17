Read("read_new.g");
Read("main_new.g");
LoadPackage( "profil" );

createMatrixFile := function( n, q, rank )
    local A, B, proj, matrixFile;
    A := RandomMat( n, n, GF(q) );
    B := Representative( GL( n, GF(q) ) );
    proj := DiagonalMat(
        One(GF(q)) * Concatenation(
            List( [1..rank], x -> 1 ),
            List( [rank+1..n], x -> 0 )
        )
    );
    A := (proj * A)^B;
    matrixFile := Concatenation( String(n), "-F", String(q), "-rk", String(rank));
    PrintTo( matrixFile, A );
    Exec( Concatenation( "mv ", matrixFile, " examples/" ) );
end;

profileGauss := function( matrixFile, chopSize  )
    local file, A, res, profile, profileName;
    profileName := Concatenation( matrixFile, "-chop", String(chopSize) );
    file := InputTextFile( Concatenation( "examples/", matrixFile ) );
    Print("Reading input...\c");
    A := EvalString( ReadAll( file ) );
    Print("OK\n");
    Exec( Concatenation( "rm -f ", profileName, ".gz" ) );
    Print("Computing...\c");
    ProfileLineByLine( Concatenation( profileName, ".gz" ) );
    res := GaussParallel( A, chopSize, chopSize );;
    UnprofileLineByLine();
    Print("OK\n");
    Print("Reading profile...\c");
    profile := ReadLineByLineProfile( Concatenation( profileName, ".gz" ) );
    Print("OK\n");
    Exec( Concatenation( "mkdir -p profiles/", profileName ) );
    Exec( Concatenation( "rm -f profiles/", profileName, "/*" ) );
    OutputAnnotatedCodeCoverageFiles( profile, Concatenation( "profiles/", profileName ) );
    return res;
end;
