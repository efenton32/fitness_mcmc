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


def lineage_count(population, bc=False):
    bc_file = get_dir(population + "d0_j.csv", "raw_data")
    bc_read = pd.read_csv(bc_file, sep=",", header=0)
    lineages = len(bc_read.index)

    if bc:
        return bc_read["Gen"].tolist()
    else:
        return lineages


def format_data(output, population, environment, replicates, bc_ass, max_days, gens_per_day):
    """
    Pipeline for formatting daily fitness assay barcode counts into a single file. Currently only works when provided a
    master file with barcode associations and some relevant naming scheme for timepoints. This will work with the first
    batch of fitness assays from Tanush.

    Parameters:
        output(str) - label for output files
        population(str) - the LTEE population, may be removed in the future
        environment(str) - growth environment of the assay
        replicates(int) - the number of total replicates to format
        bc_ass(str) - file_suffix with file type for the barcode association
        max_days(int) - number of days in the fitness assay
        gens_per_day(int) - number of generations per daily transfer of fitness assay (float support coming soon)
    """
    pop = population
    env = environment

    # Looping over replicates
    for trial in range(1, replicates + 1):
        rep = str(trial)

        header = ["BC"]
        for g in range(0, max_days + 1):
            header.append(str(g * gens_per_day))

        lineages = lineage_count(pop)

        # dict of generation value : bc count at day 0
        gens = {}

        # dict of barcode sequence : generation value
        barcodes = {}

        # final output file (column 0 is generation values, row 0 is bc counts for first ancestor over each timepoint)
        frame = np.zeros((lineages, max_days + 2))

        # list of generation values
        ind = []

        # Reading barcode association file and initializing frame
        bc_file = get_dir(population + bc_ass, "raw_data")
        with open(bc_file, 'rt') as f:
            reader = csv.reader(f)
            labels = next(reader)
            for i in range(0, lineages):
                bc = next(reader)
                barcodes[bc[0]] = int(bc[2])
                gens[int(bc[2])] = int(bc[1])
                ind.append(int(bc[2]))
        f.close()
        ind.sort()
        point = {}
        for i in range(0, len(ind)):
            frame[i][0] = ind[i]
            frame[i][1] = gens[ind[i]]
            point[ind[i]] = i

        # Reading daily data files
        for r in range(1, max_days + 1):
            list_file = pop + "_" + env + "_" + rep + "_" + str(r) + ".csv"
            data_file = get_dir(list_file, "raw_data")
            first_read = pd.read_csv(data_file, sep=",", header=0)
            lim = len(first_read.index)
            with open(data_file, "rt") as g:
                reader = csv.reader(g)
                labels = next(reader)
                for p in range(0, lim):
                    data = next(reader)
                    index = barcodes[data[0]]
                    frame[point[index]][r + 1] = int(data[1])
            g.close()

        file_out = get_dir(output + "_" + pop + "_" + env + rep + ".csv", "experimental_data")
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
