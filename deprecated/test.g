#############################################################################
##
##                             Gauss-par package
##  test.g
##                                                          Sergio Siccha
##                                                          Jendrik Brachter
##
##  Copyright...
##
##  Comparing various functions.
##  Step1:
##      - task-dependency version
##      - old sequential
##
#############################################################################

##################################################
# function compare
# Input:
#   func
#   funcAlt
#   args
#
# Output:
#   res = resAlt;
##################################################
compare := function( func, funcAlt, args )
    res := CallFuncList( func, args );
    resAlt := CallFuncList( funcAlt, args );
    return res = resAlt;
end;

randomMatrices := function( n )
    ## Create random matrices
    ## TODO
end;

testStep1 := function( iterations )
    ## Run tests
    for i in [ 1 .. iterations ] do
        equality := compare( Step1, _DEPRECATED_Step1, [ A, n ] );
        if not equality then
            ## Write A together with a hash into matrix-file
            hash := Random( 10^6, 10^7-1 );
            string := Concatenation( hash, " ", String(A) );
            AppendTo( "error-matrices.txt", string );
            Error( hash, " Task-dep and sequential Step1: different results" );
        fi;
    od;
end;

testGauss := function( iterations )
    GaussParallel(A).vectors;
    EchelonMatTransformation(A).vectors;
end;
