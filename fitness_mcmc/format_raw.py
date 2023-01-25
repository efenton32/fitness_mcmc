import csv
import numpy as np
import pandas as pd


def format_data(population, environment, replicates):
    pop = population
    env = environment

    for trial in range(1, replicates + 1):
        rep = str(trial)
        timepoints = [500, 1000, 1500, 2000, 5000, 10000, 15000, 30000, 40000, 50000]

        header = ["BC"]

        for g in range(0, 6):
            header.append(str(g * 7))

        gens = {}                       # dict of generation value:index in the frame array
        frame = np.zeros((29, 7))       # final output file (column 0 is generation values, row 0 is bc counts for first ancestor over each timepoint)
        for i in range(0, 19):
            frame[i][0] = i
            gens[i] = i
        for i in range(0, len(timepoints)):
            frame[i + 19][0] = timepoints[i]
            gens[timepoints[i]] = i + 19

        barcodes = {}                   # dict of barcode sequence: generation value
        with open("m5d0_j.csv", 'rt') as f:
            reader = csv.reader(f)
            labels = next(reader)
            for i in range(1, 30):
                bc = next(reader)
                barcodes[bc[0]] = int(bc[2])
                index = gens[barcodes[bc[0]]]
                frame[index][1] = int(bc[1])
        f.close()

        for r in range(1, 6):
            list_file = pop + "_" + env + "_" + rep + "_" + str(r) + ".csv"
            first_read = pd.read_csv(list_file, sep=",", header=0)
            lim = len(first_read.index)
            with open(list_file, "rt") as g:
                reader = csv.reader(g)
                labels = next(reader)
                for p in range(0, lim):
                    data = next(reader)
                    point = barcodes[data[0]]
                    index = gens[point]
                    frame[index][r + 1] = int(data[1])
            g.close()

        out = open("LTEE_" + pop + "_" + env + rep + ".csv", "w")
        holder = header[0]
        for w in range(1, len(header)):
            holder += "," + header[w]
        out.write(holder)
        out.write("\n")
        for l in range(0, len(frame)):
            holder = str(frame[l][0])
            for w in range(1, len(frame[l])):
                holder += "," + str(frame[l][w])
            out.write(holder)
            out.write("\n")

        out.close()
