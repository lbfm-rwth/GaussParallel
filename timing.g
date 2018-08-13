GET_REAL_TIME_OF_FUNCTION_CALL := function ( method, args, options... )
  local first_time, firstSeconds,
    firstMicroSeconds, result, second_time, secondSeconds,
  secondMicroSeconds, total, seconds, microSeconds;

  if options = [] then
    options := rec();
  else
    options := options[1];
  fi;
  if not IsBound( options.passResult ) then
    options.passResult := false;
  fi;

  first_time := IO_gettimeofday(  );
  firstSeconds := first_time.tv_sec;
  firstMicroSeconds := first_time.tv_usec;

  result := CallFuncList( method, args );

  second_time := IO_gettimeofday(  );
  secondSeconds := second_time.tv_sec;
  secondMicroSeconds := second_time.tv_usec;

  seconds := (secondSeconds - firstSeconds);
  microSeconds := secondMicroSeconds - firstMicroSeconds;
  total := seconds * 10^6 + microSeconds;
  return rec( result := result, time := total );
end;

threeSignificantDigits := function( x )
  local count;
  if not IsFloat(x) then
    Error( "x must be a Float.\n" );
  fi;
  count := 0;
  while x >= 10. do
    x := x/10;
    count := count + 1;
  od;
  # Round to three significant digits
  x := Floor( x * 100 );
  return Round(x * 10^(count-2));
end;

GetStatistics := function( data )
  local statistics;
  ## Fill the statistics vector
  Sort( data );
  # maximal string length = 12
  statistics := [
    Minimum( data ),
    0.,
    Average( data ),
    Median( data ),
    0.,
    Maximum( data )
  ];
  if Length( data ) >= 4 then
    statistics[2] := data[ Int( 0.25 * Length(data) ) ];
    statistics[5] := data[ Int( 0.75 * Length(data) ) ];
  fi;
  statistics := 1. * statistics;
  statistics := List( statistics, threeSignificantDigits );
  return statistics;
end;

## R-microbenchmark like statistics
Benchmark := function( func, args, opt... )
  local timings, columnNames, statistics, i, t, res;
  Info(InfoGauss, 1, "Start Benchmark in timing.g");
  if opt = [] then
    opt := rec();
  else
    opt := opt[1];
  fi;
  if not IsBound( opt.warmup ) then
    opt.warmup := 0;
  fi;
  if not IsBound( opt.times ) then
    opt.times := 5;
  fi;
  timings := [];
  statistics := [];

  ## Perform the computations
  if opt.warmup > 0 then
    for i in [1 .. opt.warmup] do
      GET_REAL_TIME_OF_FUNCTION_CALL( func, args );
    od;
  fi;
  for i in [ 1 .. opt.times ] do
    res := GET_REAL_TIME_OF_FUNCTION_CALL( func, args );
    Info(InfoGauss, 1, "GET_REAL_TIME_OF_FUNCTION_CALL calculation ", i);
    t := res.time;
    # We don't care about microseconds
    t := Floor( 1.0 * t / 1000 );
    timings[ Length(timings)+1 ] := t;
  od;

  statistics := GetStatistics( timings );

  res := rec( timings := timings, statistics := statistics );
  return res;
end;
