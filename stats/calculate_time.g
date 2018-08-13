# Calculate the average duration of Chief or ChiefParallel for one specific
# type of matrices and write it into the corresponding csv file.

Read("stats/type_of_matrix.g");
Read("read.g");

measuredTime := CalculateAverageTime(isParallel, dimension, dimension, rank, ring, numberChops, numberChops);

# Append obtained information to corresponding file.
if isParallel then
    AppendTo("stats/times_par.csv", dimension, ",", dimension, ",", rank, ",", ring, ",", measuredTime[3], ",", measuredTime[4], "\n");
else
    AppendTo("stats/times_seq.csv", dimension, ",", dimension, ",", rank, ",", ring, ",", measuredTime[3], ",", measuredTime[4], "\n");
fi;

FORCE_QUIT_GAP(1);
