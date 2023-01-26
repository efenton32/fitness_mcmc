# Imports
import numpy as np
import pymc3 as pm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import arviz as az
import fitness_mcmc
import fitness_mcmc.data_io as io
import fitness_mcmc.fitness_mcmc as m
import fitness_mcmc.format_raw as raw


def save_fitness(filename, lineages, s_vals, replicates, stats):
    """
    Will save the lineage labels and associated fitness values into a csv file.

    Parameters:
        filename(str) - name of output file
        lineages(list) - list of the lineage labels/barcodes
        s_vals(ndarray) - 2D numpy array (# of reps x # of lineages) of the inferred fitness values
        replicates(int) - num of replicates
        stats(ndarray) - 2D numpy aray (# of lineages x 2), contains mean and standard deviation
    """
    out_file = raw.get_dir(filename, "out")
    out = open(out_file, "w")
    header = "BC"
    for rep in range(1, replicates + 1):
        header += ",s_" + str(rep)
    out.write(header + ",average,standard deviation\n")
    for w in range(0, len(lineages)):
        line = str(lineages[w])
        for rep in range(1, replicates + 1):
            line += "," + str(s_vals[rep - 1][w])
        line += "," + str(stats[w][0]) + "," + str(stats[w][1])
        out.write(line + "\n")
    out.close()


def fitness_pipeline(population, environment, reference=1, replicates=3, format=False):
    """
    Serves as the master pipeline for getting a single environment's (and replicates) inferred fitness values. Will
    update these comments as I develop functionality here. Currently will:
    -Will handle multiple replicates
    -Format raw barcode count files to be compatible with data_io & fitness_mcmc
    -Calculate fitness inferences from fitness_mcmc via find_MAP
    -Save the original versus reconstructed frequency plots for each replicate
    -Save fitness values for all replicates into a single csv file with mean and standard deviation.

    Parameters:
        population(str) - LTEE population of interest
        environment(str) - environment of fitness assay
        reference(int) - reference lineage defined to have s = this value
        replicates(int) - # of replicates to be processed
        format(bool) - whether or not to format from raw barcode counts
    """
    if format:
        raw.format_data(population, environment, replicates)

    # Loop for all replicates
    name = "LTEE_" + population + "_" + environment
    barcodes = []
    lin_count = raw.lineage_count("pop1")
    s = np.zeros((replicates, lin_count))
    for rep in range(1, replicates + 1):
        s[rep - 1][0] = reference
        # Loading data with data_io
        r = name + str(rep)
        data, time, ordered_counts = io.load_data(r + ".csv", return_ordered=False, delimiter=",")
        if rep == 1:
            barcodes = data["BC"].tolist()

        # Creating fitness model and getting inferences with fitness_mcmc
        fitness_model = m.Fitness_Model(ordered_counts, time, s_ref=reference, prior="flat")
        fitness_model.find_MAP()

        # Retrieving "s" and "f0" values from the model and recording
        raw_s = fitness_model.map_estimate["s"]
        for lin in range(0, len(raw_s)):
            s[rep - 1][lin + 1] = raw_s[lin][0]

        # Plotting MAP
        output = raw.get_dir(r + "_freq.png", "out")
        fitness_model.plot_MAP_estimate(type="lin", filename=output)

    # Getting replicate average and standard deviation
    stats = np.zeros((lin_count, 2))
    for x in range(0, lin_count):
        reps = []
        for q in range(0, replicates):
            reps.append(s[q][x])
        stats[x][0] = np.mean(reps)
        stats[x][1] = np.std(reps)

    # Writing output
    save_fitness(name + "_fitness.csv", barcodes, s, replicates, stats)
    # Will add plotting function here
