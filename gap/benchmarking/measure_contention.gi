GAUSS_prepareCounters := function()
    Sleep(1);
    THREAD_COUNTERS_ENABLE();
    THREAD_COUNTERS_RESET();
    return ThreadID(CurrentThread());
end;

GAUSS_getCountersForTaskThreads := function()
    local counters;
    Sleep(1);
    counters := [ThreadID(CurrentThread())];
    Append(counters, THREAD_COUNTERS_GET());
    return counters;
end;

GAUSS_prepare := function(nrAvailableThreads)
    local tasks, taskNumbers;
    tasks := List([1..nrAvailableThreads],
                  x -> RunTask(GAUSS_prepareCounters));
    taskNumbers := List(tasks, TaskResult);
    if Size(Set(taskNumbers)) <> GAPInfo.KernelInfo.NUM_CPUS then
        ErrorNoReturn("GaussPar: Unable to activate counters for all threads");
    fi;
end;

GAUSS_compute := function(A, q, numberBlocks, showOutput)
    local res;
    if (showOutput) then
        Print("Starting the parallel algorithm:\n",
              "    EchelonMatTransformationBlockwise.\n\n");
    fi;
    res := GAUSS_GET_REAL_TIME_OF_FUNCTION_CALL(
        EchelonMatTransformationBlockwise,
        [A, rec( galoisField := GF(q), 
                 numberBlocks := numberBlocks )]);
    # The variable `time` stores the cpu time that was used by the execution
    # of the last statement.
    res.cpuTime := time;
    return res;
end;

GAUSS_evaluate := function(nrAvailableThreads, bench, A, showOutput)
    local CPUTimeCompute, resPar, benchStd, resStd, correct, counters,
        counter, totalAcquired, totalContended, ratio;

    resPar := bench.result;
    # Use ms instead of microseconds
    bench.time := Round(bench.time / 1000.);
    if (showOutput) then
        Print("Wall time parallel execution (ms): ", bench.time, "\n");
    fi;
    counters := List([1..nrAvailableThreads],
                     x -> RunTask(GAUSS_getCountersForTaskThreads));
    counters := List(counters, TaskResult);
    SortParallel(List(counters, x -> x[1]), counters);
    totalAcquired := Sum(List(counters, x -> x[2]));
    totalContended := Sum(List(counters, x -> x[3]));
    ratio := Round(totalContended / totalAcquired * 100.);
    if (showOutput) then
        Print("Lock statistics(estimates):\n");
        Print("acquired - ", totalAcquired, ", contended - ", totalContended);
        Print(", ratio - ", ratio, "%\n");
        Print("Locks acquired and contention counters per thread\n");
        Print("[ thread, locks acquired, locks contended ]:\n");
        for counter in counters do
            Print(counter, "\n");
        od;
        Print("\n");
        Print("Starting the sequential algorithm:\n");
        Print("    EchelonMatTransformation.\n\n");
    fi;
    benchStd := GAUSS_GET_REAL_TIME_OF_FUNCTION_CALL(EchelonMatTransformation, [A]);
    # Use ms instead of microseconds
    benchStd.time := Round(benchStd.time / 1000.);
    resStd := benchStd.result;
    correct := resStd.vectors = resPar.vectors
           and resStd.coeffs = resPar.coeffs;
    if not correct then
        ErrorNoReturn("GaussPar: Result incorrect!");
    fi;
    if (showOutput) then
        Print("Wall time Gauss pkg execution (ms): ", benchStd.time, "\n\n");
        Print("Speedup factor (sequential / parallel wall time):\n");
        Print(GAUSS_threeSignificantDigits(1. * benchStd.time / bench.time));
        Print("\n");
    fi;
end;

InstallGlobalFunction("GAUSS_MeasureContention",
function(numberBlocks, q, A, options...)
    local nrAvailableThreads, bench, showOutput;

    # showOutput can be set for testing purposes
    showOutput := true;
    if ((Length(options) = 1) and (options[1] = false)) then
        showOutput := false;
    fi;
    if showOutput then
        Print("Make sure you called GAP with sufficient preallocated",
              " memory via `-m` if you try non-trivial examples!",
              " Otherwise garbage collection will be a big overhead.\n\n");
    fi;

    nrAvailableThreads := GAPInfo.KernelInfo.NUM_CPUS;
    bench := "";

    GAUSS_prepare(nrAvailableThreads);;
    bench := GAUSS_compute(A, q, numberBlocks, showOutput);;
    GAUSS_evaluate(nrAvailableThreads, bench, A, showOutput);;
end);
