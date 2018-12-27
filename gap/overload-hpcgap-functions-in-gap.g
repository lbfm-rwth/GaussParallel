# Mock-ups for tasks in GAP
ScheduleTask := function(cond, func, args...)
    return CallFuncListWrap(RunTask, Concatenation([func], args))[1];
end;

WaitTask := function(args...)
    return;
end;

BindGlobal("RegionOf", x -> 0);
