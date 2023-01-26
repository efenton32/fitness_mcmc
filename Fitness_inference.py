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


def fitness_pipeline(population, environment, output, replicates=3, format=True):
    """
    Serves as the master pipeline for getting a single environment's (and replicates) inferred fitness values. Will
    update these comments as I develop functionality here.

    Parameters:
        population(str) - LTEE population of interest
        environment(str) - environment of fitness assay
        output(str) - tag for output files
        replicates(int) - # of replicates to be processed
        format(bool) - whether or not to format from raw barcode counts
    """
    if format:
        raw.format_data(population, environment, replicates)

    # Loop for all replicates
    for rep in range(1, replicates + 1):
        # Loading data with data_io
        r = population + "_" + environment + str(rep)
        data, time, ordered_counts = io.load_data("LTEE_" + r + ".csv", return_ordered=False, delimiter=",")

        # Creating fitness model and getting inferences with fitness_mcmc
        fitness_model = m.Fitness_Model(ordered_counts, time, s_ref=1, prior="flat")
        fitness_model.find_MAP()

        # Retrieving "s" and "f0" values from the model
        raw_s = fitness_model.map_estimate["s"]
        vals_f0 = fitness_model.map_estimate["f0"]
        vals_s = np.zeros((29))
        vals_s[0] = 0
        for i in range(0, len(raw_s)):
            vals_s[i + 1] = raw_s[i]
        vals_data = fitness_model.data


