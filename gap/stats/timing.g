## Calculates time statistics for one matrix of a specific type using Benchmark()
GAUSS_CalculateTime := function(isParallel, height, width, rank, ring, numberChopsH, numberChopsW, randomSeed)
    local echelon, shapeless, result, times, r;
    Info(InfoGauss, 1, "Start CalculateTime");
    times := 0;

    # Create random matrices, calculate time.
    echelon := RandomEchelonMat(height, width, rank, randomSeed, ring);;
    Info(InfoGauss, 3, "Echelon matrix:");
    Info(InfoGauss, 3, echelon);
    shapeless := GAUSS_shapelessMat(echelon, height, width, randomSeed, ring);;
    Info(InfoGauss, 3, "Shapeless matrx:");
    Info(InfoGauss, 3, shapeless);
    if isParallel then
        Info(InfoGauss, 2, "Parallel version:");
    else
        Info(InfoGauss, 2, "Sequential version:");
    fi;
    times := GAUSS_Benchmark(DoEchelonMatTransformationBlockwise, [shapeless, ring, isParallel, numberChopsH, numberChopsW]);

    return times.timings;
end;

## Calculates time statistics for 10 matrices of a specific type
GAUSS_CalculateAverageTime := function(isParallel, height, width, rank, ring, numberChopsH, numberChopsW)
    local randomSeed, timings, statistics, i;

    Info(InfoGauss, 1, "Start CalculateAverageTime in stats/timing.g");

    randomSeed := RandomSource(IsMersenneTwister);;
    timings := [];

    # Do a few times and calculate average.
    for i in [ 1 .. 10 ] do
        Info(InfoGauss, 2, "CalculateTime calculation no.", i);
        Append(timings, GAUSS_CalculateTime(isParallel, height, width, rank, ring, numberChopsH, numberChopsW, randomSeed));
    od;

    statistics := GAUSS_GetStatistics(timings);
    return statistics;
end;
