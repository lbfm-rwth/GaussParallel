Read("read_new.g");
Read("main_new.g");
LoadPackage( "profil" );

profileGauss := function( n, q, chopSize  )
    local A, res, profile;
    A := RandomMat( n, n, GF(q) );
    Exec("rm profile.gz");
    ProfileLineByLine( "profile.gz" );
    res := GaussParallel( A, chopSize, chopSize );;
    UnprofileLineByLine();
    profile := ReadLineByLineProfile( "profile.gz" );
    Exec("mkdir profile");
    Exec("rm -r profile/*");
    OutputAnnotatedCodeCoverageFiles( profile, "profile" );
    return res;
end;
