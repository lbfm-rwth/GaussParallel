Read("read_hpc.g");

#Make sure you called GAP with sufficient
#preallocated memory via `-m`

prepareCounters := function()
    Sleep(1);
    THREAD_COUNTERS_ENABLE();
    THREAD_COUNTERS_RESET();
    return ThreadID(CurrentThread());
end;

getCountersForTaskThreads := function()
    local counters;
    Sleep(1);
    counters := [ThreadID(CurrentThread())];
    Append(counters, THREAD_COUNTERS_GET());
    return counters;
end;

prepare := function(nrAvailableThreads)
    local tasks, taskNumbers;
    tasks := List([1..nrAvailableThreads],
                  x -> RunTask(prepareCounters));
    taskNumbers := List(tasks, TaskResult);
    if Size(Set(taskNumbers)) <> GAPInfo.KernelInfo.NUM_CPUS then
        ErrorNoReturn("GaussPar: Unable to activate counters for all threads");
    fi;
end;

compute := function(A, q, numberChops)
    return GET_REAL_TIME_OF_FUNCTION_CALL(ChiefParallel,
                                            [GF(q), A, numberChops, numberChops]);
end;

evaluate := function(nrAvailableThreads, bench, A, visible)
    local CPUTimeCompute, resPar, benchStd, resStd, correct, counters,
    totalAcquired, totalContended, factor;

    CPUTimeCompute := time;
    resPar := bench.result;
    benchStd := GET_REAL_TIME_OF_FUNCTION_CALL(EchelonMatTransformation, [A]);
    resStd := benchStd.result;
    correct := -1 * resStd.vectors = resPar.vectors
           and -1 * resStd.coeffs = resPar.coeffs;
    if not correct then
        ErrorNoReturn("GaussPar: Result incorrect!");
    fi;
    counters := List([1..nrAvailableThreads], x -> RunTask(getCountersForTaskThreads));
    counters := List(counters, TaskResult);
    SortParallel(List(counters, x -> x[1]), counters);
    totalAcquired := Sum(List(counters, x -> x[2]));
    totalContended := Sum(List(counters, x -> x[3]));
    factor := Round(totalContended / totalAcquired * 100.);
    
    if (visible) then
        Print("Wall time  parallel execution: ", bench.time, "\n");
        Print("CPU  time  parallel execution: ", CPUTimeCompute*1000, "\n");
        Print("Wall time Gauss pkg execution: ", benchStd.time, "\n");
        Print("Lock statistics(estimates):\n");
        Print("acquired - ", totalAcquired, ", contended - ", totalContended);
        Print(", factor - ", factor, "%\n");
        Print("Locks acquired and contention counters per thread:\n");
        Print(counters, "\n");
    fi;
end;

MeasureContention := function(numberChops, q, A, options...)
    local nrAvailableThreads, bench, visible;

    visible := true;
    if ((Length(options) = 1) and (options[1] = false)) then
        visible := false;
    fi;

    nrAvailableThreads := GAPInfo.KernelInfo.NUM_CPUS;
    bench := "";
    
    prepare(nrAvailableThreads);;
    bench := compute(A, q, numberChops);;
    evaluate(nrAvailableThreads, bench, A, visible);;
end;
