#!/usr/bin/env python3
# Creates two csv files and writes the statistics for the sequential and
# parallel version of the gaussian algorithm into those files.
import subprocess
import itertools
import sys
import argparse

# TODO: --small suite flag
# round numbers

parser = argparse.ArgumentParser(
    description = "gap and hpcgap have to be available in the path or paths to them must be passed as arguments"
)
parser.add_argument("--debug",
                    help = "only show one set of commands and don't pass it to gap/hpcgap",
                    action = "store_true")
parser.add_argument("--small-suite",
                    help = "only run a small subset of the suite",
                    action = "store_true")
command_line_options = parser.parse_args()


# Write files.
subprocess.call('echo "height,width,rank,ring,average,median" > stats/times_par.csv', shell=True)
subprocess.call('echo "height,width,rank,ring,time,average,median" > stats/times_seq.csv', shell=True)

# Calculates average and median of the duration of Chief for every combination
# from lists of width, height, rank, ring. Saves them in the csv files.
# TODO: make ranks and numberBlocks depend on the dimension
if not command_line_options.small_suite:
    isParallel = ["true", "false"]
    dimensions = [5, 10, 50, 100, 200, 500]
    rings = [2, 3, 5, 11, 17]
    ranks = [1, 3, 7, 15, 50, 90, 155, 350, 500]
    numberBlocks = [1, 5, 50]
else:
    isParallel = ["true", "false"]
    dimensions = [5, 10]
    rings = [2, 3, 5]
    ranks = [1, 3, 5, 7]
    numberBlocks = [1, 2, 3]

specifications = list(itertools.product(isParallel, dimensions, ranks, rings, numberBlocks))
for (p, d, rank, ring, numberBlocks) in specifications:
    if rank <= d and numberBlocks < d:
        print(p, d, rank, ring, numberBlocks)
        if p == "true":
            outfile = '"stats/times_par.csv"'
            gap = 'hpcgap'
        else:
            outfile = '"stats/times_seq.csv"'
            gap = 'gap-master'
        args = 'isParallel := ' + p \
            + ';; dimension :=' + str(d) \
            + ';; rank := ' + str(rank) \
            + ';; ring := GF(' + str(ring) + ')' \
            + ';; numberBlocks := ' + str(numberBlocks) \
            + ';;\n'
        instructions = "\n".join([
            'LoadPackage("GaussPar");;',
            'LoadPackage("io");;',
            'ReadPackage("GaussPar", "gap/benchmarking/timing.g");;',
            'measuredTime := GAUSS_CalculateAverageTime(isParallel, '
            + 'dimension, dimension, rank, ring, numberBlocks, numberBlocks);;',
            'AppendTo(' + outfile 
            + ', dimension, ",", dimension, ",", rank, ",", ring, ",", '
            + 'measuredTime[3], ",", measuredTime[4], "\\n");'
        ])
        if command_line_options.debug:
            print(args)
            print(instructions)
            sys.exit(0)
        subprocess.call(gap + ' << EOF\n'
            + args
            + instructions
            + '\nEOF\n'
            , shell=True)
