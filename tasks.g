# Mock-ups for tasks in GAP
ScheduleTask := function(cond, func, args...)
    return RunTask(func, args);
end;

WaitTask := function(args...)
    return;
end;
