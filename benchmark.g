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

  if options.passResult then
    result := CallFuncList( method, args );
  else
    CallFuncList( method, args );
  fi;

  second_time := IO_gettimeofday(  );
  secondSeconds := second_time.tv_sec;
  secondMicroSeconds := second_time.tv_usec;

  seconds := (secondSeconds - firstSeconds);
  microSeconds := secondMicroSeconds - firstMicroSeconds;
  total := seconds * 10^6 + microSeconds;
  if options.passResult then
    return rec( result := result, time := total );
  else
    return total;
  fi;
end;

alignRight := function( obj, n )
  ## fill up string with " " on the left side
  local m, string;
  string := String( obj );
  m := Maximum( 0, n - Size(string) );
  return Concatenation(
    RepeatedString( " ", m ),
    string
  );
end;

alignCenter := function( obj, n )
  ## fill up string with " " on both sides
  local m, ml, mr, string;
  string := String( obj );
  m := Maximum( 0, n - Size(string) );
  if IsEvenInt( m ) then
    ml := m/2;
    mr := m/2;
  else
    ml := Int( Floor( 1. * m/2 ) );
    mr := Int( Floor( 1. * m/2 ) ) + 1;
  fi;
  return Concatenation(
    RepeatedString( " ", ml ),
    string,
    RepeatedString( " ", mr )
  );
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
  # Rat: Cast to Rationals to avoid FLOPS errors
  x := Floor( x * 100 );
  return x * 10^(count-2);
end;

## R-microbenchmark like Printing of statistics
PrintStats := function( statistics )
  ## TODO enable multiple row statistics and printing of these
  ## TODO first determine timeUnit, then round to three significant digits appropriately
  local timeUnit, median, columnNames, count;
  statistics := List( statistics, threeSignificantDigits );
  columnNames := [ "Min", "1st Quart", "Mean", "Median", "3rd Quart", "Max" ];
  median := Median( statistics );
  if median < 10^4 * 1. then
    timeUnit := "microseconds";
  elif median < 10^7 * 1. then
    timeUnit := "milliseconds";
    statistics := statistics / (10^3);
  elif median < 10^10 * 1. then
    timeUnit := "seconds";
    statistics := statistics / (10^6);
  elif median < 60 * 10^10 * 1. then
    timeUnit := "minutes";
    statistics := statistics / (60 * 10^6);
  else
    timeUnit := "hours";
    statistics := statistics / (60 * 60 * 10^6);
  fi;
  Print( timeUnit, "\n" );
  columnNames := List( columnNames, x -> alignRight( x, 13 ) );
  columnNames := Concatenation( columnNames );
  Print( columnNames, "\n" );
  statistics := List( statistics, x -> alignRight( Int(x), 13 ) );
  statistics := Concatenation( statistics );
  Print( statistics, "\n" );
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
  return 1. * statistics;
end;

## R-microbenchmark like statistics
Benchmark := function( method, args, opt... )
  local timings, columnNames, statistics, i, t, res;
  if opt = [] then
    opt := rec();
  else
    opt := opt[1];
  fi;
  if not IsBound( opt.warmup ) then
    opt.warmup := 5;
  fi;
  if not IsBound( opt.times ) then
    opt.times := 100;
  fi;
  if not IsBound( opt.withResult ) then
    opt.withResult := false;
  fi;
  timings := [];
  statistics := [];

  ## Perform the computations
  if opt.warmup > 0 then
    for i in [1 .. opt.warmup] do
      t := GET_REAL_TIME_OF_FUNCTION_CALL( method, args );
    od;
  fi;
  for i in [ 1 .. opt.times ] do
    t := GET_REAL_TIME_OF_FUNCTION_CALL( method, args );
    timings[ Length(timings)+1 ] := 1.0 * t;
  od;

  statistics := GetStatistics( timings );

  PrintStats( statistics );
  res := rec( timings := timings, statistics := statistics );
  if opt.withResult then
    res.res := CallFuncList( method, args );
  fi;
  return res;
end;

## Execute function call only once
BenchmarkOnce := function( method, args, opt... )
  return Benchmark( method, args, rec( warmup := 0, times := 1  ) );
end;
