import os
import csv
import numpy as np
import pandas as pd


def get_dir(file_name, sub_dir):
    """
    This function will return the correct path to the relevant data files for the LTEE data.
    This file is designed to work with the processed barcode counts from LTEE fitness assays. For other users,
    please format the data as seen in the regular documentation.

    Parameters:
        file_name(str) - name of the data file including file type
        sub_dir(str) - name of the folder in which the desired file resides

    Returns:
        data_path(str) - the exact path to where the named file resides
    """
    start = os.path.abspath(__file__)
    start_dir = os.path.dirname(start)
    data_path = os.path.join(start_dir, sub_dir, file_name)

    return data_path


def format_data(population, environment, replicates):
    """
    Pipeline for formatting daily fitness assay barcode counts into a single file. Currently only works when provided a
    master file with barcode associations and some relevant naming scheme for timepoints. This will work with the first
    batch of fitness assays from Tanush, but will need to be updated for robustness in the future. For now, we ball.

    Parameters:
        population(str) - the LTEE population, may be removed in the future
        environment(str) - growth environment of the assay
        replicates(int) - the number of total replicates to format
    """
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
        bc_file = get_dir("m5d0_j.csv", "raw_data")
        with open(bc_file, 'rt') as f:
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
            data_file = get_dir(list_file, "raw_data")
            first_read = pd.read_csv(data_file, sep=",", header=0)
            lim = len(first_read.index)
            with open(data_file, "rt") as g:
                reader = csv.reader(g)
                labels = next(reader)
                for p in range(0, lim):
                    data = next(reader)
                    point = barcodes[data[0]]
                    index = gens[point]
                    frame[index][r + 1] = int(data[1])
            g.close()

        file_out = get_dir("LTEE_" + pop + "_" + env + rep + ".csv", "experimental_data")
        out = open(file_out, "w")
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
        print("formatting done")
