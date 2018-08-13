## Calculates time statistics for one matrix of a specific type using Benchmark()
CalculateTime := function(isParallel, height, width, rank, ring, numberChopsH, numberChopsW, randomSeed)
    local echelon, shapeless, result, times, r;
    Info(InfoGauss, 1, "Start CalculateTime");
    times := 0;

    # Create random matrices, calculate time.
    echelon := RandomEchelonMat(height, width, rank, randomSeed, ring);;
    Info(InfoGauss, 3, "Echelon matrix:");
    Info(InfoGauss, 3, echelon);
    shapeless := _GAUSS_shapelessMat(echelon, height, width, randomSeed, ring);;
    Info(InfoGauss, 3, "Shapeless matrx:");
    Info(InfoGauss, 3, shapeless);
    if isParallel then
        Info(InfoGauss, 2, "Parallel version:");
        times := Benchmark(ChiefParallel, [ring, shapeless, numberChopsH, numberChopsW]);
    else
        Info(InfoGauss, 2, "Sequential version:");
        times := Benchmark(Chief, [ring, shapeless, numberChopsH, numberChopsW]);
    fi;

    return times.timings;
end;

## Calculates time statistics for 10 matrices of a specific type
CalculateAverageTime := function(isParallel, height, width, rank, ring, numberChopsH, numberChopsW)
    local randomSeed, timings, statistics, i;

    Info(InfoGauss, 1, "Start CalculateAverageTime in stats/timing.g");

    randomSeed := RandomSource(IsMersenneTwister);;
    timings := [];

    # Do a few times and calculate average.
    for i in [ 1 .. 10 ] do
        Info(InfoGauss, 2, "CalculateTime calculation no.", i);
        Append(timings, CalculateTime(isParallel, height, width, rank, ring, numberChopsH, numberChopsW, randomSeed));
    od;

    statistics := GetStatistics(timings);
    return statistics;
end;
