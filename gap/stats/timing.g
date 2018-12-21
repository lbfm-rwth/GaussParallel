## Calculates time statistics for one matrix of a specific type using Benchmark()
GAUSS_CalculateTime := function(isParallel, height, width, rank, ring, numberBlocksH, numberBlocksW, randomSeed)
    local echelon, shapeless, result, times, r;
    Info(InfoGauss, 2, "Start CalculateTime");
    times := 0;

    # Create random matrices, calculate time.
    echelon := RandomEchelonMat(height, width, rank, randomSeed, ring);;
    Info(InfoGauss, 4, "Echelon matrix:");
    Info(InfoGauss, 4, echelon);
    shapeless := GAUSS_RandomMatFromEchelonForm(echelon, height, width, randomSeed, ring);;
    Info(InfoGauss, 4, "Shapeless matrx:");
    Info(InfoGauss, 4, shapeless);
    if isParallel then
        Info(InfoGauss, 3, "Parallel version:");
    else
        Info(InfoGauss, 3, "Sequential version:");
    fi;
    times := GAUSS_Benchmark(
        DoEchelonMatTransformationBlockwise,
        [
            shapeless,
            rec( galoisField := ring, IsHPC := isParallel,
            numberBlocksHeight := numberBlocksH,
            numberBlocksWidth := numberBlocksW )
        ]
    );

    return times.timings;
end;

## Calculates time statistics for 10 matrices of a specific type
GAUSS_CalculateAverageTime := function(isParallel, height, width, rank, ring, numberBlocksH, numberBlocksW)
    local randomSeed, timings, statistics, i;

    Info(InfoGauss, 2, "Start CalculateAverageTime in stats/timing.g");

    randomSeed := RandomSource(IsMersenneTwister);;
    timings := [];

    # Do a few times and calculate average.
    for i in [ 1 .. 10 ] do
        Info(InfoGauss, 3, "CalculateTime calculation no.", i);
        Append(timings, GAUSS_CalculateTime(isParallel, height, width, rank, ring, numberBlocksH, numberBlocksW, randomSeed));
    od;

    statistics := GAUSS_GetStatistics(timings);
    return statistics;
end;
