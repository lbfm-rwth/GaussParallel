Read("read_hpc.g");

ExampleGaussTrafo := function(n,f)
    local A,tmp,test;
    
    Print("Generating test matrix","\n"); 
    A := RandomMat( n,n,f );
    Print("Running EchelonMatTransformation...","\n"); 
    tmp := EchelonMatTransformation(A);
    Print("Running a version of GaussTrafo...","\n"); 
    test := GaussTrafoSemiHPC( A,10,10,f ); ## Matrix,chopsizes,groundfield
    #test := GaussTrafoHPC( A,10,10,f ); ## Matrix,chopsizes,groundfield

    Print( "Equal remnant: " ,(tmp.vectors = -test.vectors) , "\n" );
    Print( "Equal transformation: " , (tmp.coeffs = -test.coeffs) , "\n" );
    
    return 0;
end;
