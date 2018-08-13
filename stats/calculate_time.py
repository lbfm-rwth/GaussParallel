# Creates two csv files and writes the statistics for the sequential and
# parallel version into those files.
import subprocess
import itertools

# Write files.
subprocess.call('echo "height,width,rank,ring,average,median" > stats/times_par.csv', shell=True)
subprocess.call('echo "height,width,rank,ring,time,average,median" > stats/times_seq.csv', shell=True)

# Calculates average and median of the duration of Chief and ChiefParallel for every combination
# from lists of width, height, rank, ring. Save them in the csv files.
# Right now the gap program that calculates the average times gets the
# specifications through reading a file that we create here. That is thread-
# safe because only ChiefParallel uses parallel programming and this file
# here doesn't. But obviously it is not nice. Change if you know something better!
isParallel = ["true", "false"]
dimensions = [5, 10, 50, 100, 200, 500]
ranks = [1, 3, 7, 15, 50, 90, 155, 350, 500]
rings = [2, 3, 5, 11, 17]
numberChops = [1, 5, 50]

specifications = list(itertools.product(isParallel, dimensions, ranks, rings, numberChops))
specifications = specifications[0:10]
for (p, d, ra, ri, n) in specifications:
    print(p, d, ra, ri, n)
    if ra <= d and n < d:
        subprocess.call('echo "isParallel := ' + p
            + ';; dimension :=' + str(d)
            + ';; rank := ' + str(ra)
            + ';; ring := GF(' + str(ri)
            + ');; numberChops := ' + str(n)
            + ';;" > stats/type_of_matrix.g',
            shell=True)
        subprocess.call('../../hpcgap/bin/gap.sh stats/calculate_time.g', shell=True) # CHANGE PATH HERE!
