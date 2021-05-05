if IsHPCGAP then
    # Below we define mock-ups for tasks in GAP. To not define e.g. a
    # ScheduleTask in GAP we use the name GAUSS_ScheduleTask both in GAP and in
    # HPC-GAP.
    GAUSS_ScheduleTask := ScheduleTask;
    GAUSS_WaitTask := WaitTask;
else
    GAUSS_ScheduleTask := function(cond, func, args...)
        return CallFuncListWrap(RunTask, Concatenation([func], args))[1];
    end;

    GAUSS_WaitTask := function(args...)
        return;
    end;

    # Add this, so that it is available for tab-completion
    GAUSS_MeasureContention := x-> "Only available in GAP";
fi;
