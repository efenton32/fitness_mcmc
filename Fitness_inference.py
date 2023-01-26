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


def save_fitness(filename, lineages, s_vals, reference):
    """
    Will save the lineage labels and associated fitness values into a csv file.

    Parameters:
        filename(str) - name of output file
        lineages(list) - list of the lineage labels/barcodes
        s_vals(array) - 2D numpy array (# of reps x # of lineages) of the inferred fitness values, excluding the reference
        reference(int) - reference fitness value
    """
    out_file = raw.get_dir(filename, "out")
    out = open(out_file, "w")
    out.write("BC,s\n")
    out.write(str(lineages[0]) + "," + str(reference) + "\n")
    for w in range(0, len(s_vals)):
        out.write(str(lineages[w + 1]) + "," + str(s_vals[w][0]) + "\n")
    out.close()


def fitness_pipeline(population, environment, reference=1, replicates=3, format=True):
    """
    Serves as the master pipeline for getting a single environment's (and replicates) inferred fitness values. Will
    update these comments as I develop functionality here.

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

    for rep in range(1, replicates + 1):
        # Loading data with data_io
        r = name + str(rep)
        data, time, ordered_counts = io.load_data(r + ".csv", return_ordered=False, delimiter=",")
        if rep == 1:
            barcodes = data["BC"].tolist()

        # Creating fitness model and getting inferences with fitness_mcmc
        fitness_model = m.Fitness_Model(ordered_counts, time, s_ref=reference, prior="flat")
        fitness_model.find_MAP()

        # Retrieving "s" and "f0" values from the model and recording
        vals_s = fitness_model.map_estimate["s"]

        # Plotting MAP
        output = raw.get_dir(r + "_freq.png", "out")
        fitness_model.plot_MAP_estimate(type="lin", filename=output)

    save_fitness(name + "_fitness.csv", barcodes, , reference)



