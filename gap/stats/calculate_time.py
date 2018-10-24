# Creates two csv files and writes the statistics for the sequential and
# parallel version of the gaussian algorithm into those files.
import subprocess
import itertools
import sys

# If we only want to check whether the script works, we dont run the algorithm
# for all possible input combinations.
if len(sys.argv) == 2:
    debugging = sys.argv[1]
    if debugging == "--debug":
        debugging = True
    else:
        print 'Available command line options: --debug'
        sys.exit(1)

# Write files.
subprocess.call('echo "height,width,rank,ring,average,median" > stats/times_par.csv', shell=True)
subprocess.call('echo "height,width,rank,ring,time,average,median" > stats/times_seq.csv', shell=True)

# Path to hpcgap executable
hpcgap = '/home/sergio/projects/gap-master/build/hpcgap/bin/gap.sh'

# Calculates average and median of the duration of Chief for every combination
# from lists of width, height, rank, ring. Saves them in the csv files.
isParallel = ["true", "false"]
dimensions = [5, 10, 50, 100, 200, 500]
ranks = [1, 3, 7, 15, 50, 90, 155, 350, 500]
rings = [2, 3, 5, 11, 17]
numberChops = [1, 5, 50]

specifications = list(itertools.product(isParallel, dimensions, ranks, rings, numberChops))
if debugging:
    specifications = specifications[0:10]
for (p, d, ra, ri, n) in specifications:
    print(p, d, ra, ri, n)
    if ra <= d and n < d:
        if p == "true":
            outfile = '"stats/times_par.csv"'
        else:
            outfile = '"stats/times_seq.csv"'
        args = 'isParallel := ' + p \
            + ';; dimension :=' + str(d) \
            + ';; rank := ' + str(ra) \
            + ';; ring := GF(' + str(ri) + ')' \
            + ';; numberChops := ' + str(n) \
            + ';;\n'
        instructions = 'measuredTime := _GAUSS_CalculateAverageTime(isParallel, ' \
            + 'dimension, dimension, rank, ring, numberChops, numberChops);;' \
            + '\n' \
            + 'AppendTo(' + outfile \
            + ', dimension, ",", dimension, ",", rank, ",", ring, ",", ' \
            + 'measuredTime[3], ",", measuredTime[4], "\\n");' \
            + '\n'
        subprocess.call(hpcgap + ' -q read.g << EOF\n'
            + args
            + instructions
            + 'EOF\n'
            , shell=True)
