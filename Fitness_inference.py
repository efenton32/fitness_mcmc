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


def read_s(out, population, environment, replicates=3):
    """
    Reads fitness_pipeline fitness files from out folder back into data (easier/faster than re-running the pipeline)

    Parameters:
        out(str) - file tag user imputs
        population(str) - population of interest
        environment(str) - fitness assay environment
    """
    name = out + "_" + population + "_" + environment + "_fitness.csv"
    data_file = raw.get_dir(name, "out")
    bc_read = pd.read_csv(data_file, sep=",", header=0)
    lineages = len(bc_read.index)
    s = np.zeros((replicates, lineages))
    stats = np.zeros((lineages, 2))
    barcodes = bc_read["BC"].tolist()
    for rep in range(0, replicates):
        k = "s_" + str(rep + 1)
        s[rep] = bc_read[k]
    stats[:, 0] = bc_read["average"]
    stats[:, 1] = bc_read["standard deviation"]

    return s, stats, barcodes


def fitness_plot(gens, s, stats, filename):
    """
    Plots LTEE lineage fitness per generation.

    Parameters:
        gens(list) - lineage labels, must be ints and distinguish ancestors
        s(array) - 2D numpy fitness array
        stats(array) - 2D numpy stats array
        filename(str) - output file name
    """
    vals = stats[:, 0]
    errs = stats[:, 1]
    plt.errorbar(gens, vals, errs, fmt='k:', capsize=3.0)
    for i in range(0, 3):
        plt.scatter(gens, s[i])
    plt.title('Inferred Fitness')
    plt.xlabel('Time (generations)')
    plt.ylabel('Fitness')
    out_file = raw.get_dir(filename, "out")
    plt.savefig(out_file)
    plt.close()


def fitness_pipeline(output, population, environment, reference=1, replicates=3, max_days=5, gens_per_day=7,
                     format=False, doc=False, recon=False):
    """
    Serves as the master pipeline for getting a single environment's (and replicates) inferred fitness values. Will
    update these comments as I develop functionality here. Currently will:
    -Will handle multiple replicates
    -Format raw barcode count files to be compatible with data_io & fitness_mcmc
    -Calculate fitness inferences from fitness_mcmc via find_MAP
    -Save the original versus reconstructed frequency plots for each replicate
    -Save fitness values for all replicates into a single csv file with mean and standard deviation.

    Parameters:
        output(str) - label for output files
        population(str) - LTEE population of interest
        environment(str) - environment of fitness assay
        reference(int) - reference lineage defined to have s = this value
        replicates(int) - # of replicates to be processed
        max_days(int) - # of days in fitness assay
        gens_per_day(int) - # of generations per day during assay (float support coming soon)
        format(bool) - whether or not to format from raw barcode counts
        doc(bool) - whether or not this is being run for documentation in a jupyter notebook
    """
    name = output + "_" + population + "_" + environment
    if format:
        raw.format_data(output, population, environment, replicates, max_days, gens_per_day)

    # Loop for all replicates
    barcodes = []
    lin_count = raw.lineage_count(name + "_1.csv")
    s = np.zeros((replicates, lin_count))
    for rep in range(1, replicates + 1):
        s[rep - 1][0] = reference
        # Loading data with data_io
        r = name + "_" + str(rep)
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
        if doc:
            fitness_model.plot_MAP_estimate(type="lin")
        elif recon:
            new_traj = fitness_model.plot_MAP_estimate(type="lin", recon=recon)
            out_traj_name = raw.get_dir(r + "_reconstructed_freq.csv", "out")
            out_traj = open(out_traj_name, "w")
            header = "BC"
            for t in range(0, max_days + 1):
                header += "," + str(t * gens_per_day)
            out_traj.write(header + "\n")
            for bc in range(0, len(barcodes)):
                holder = str(barcodes[bc])
                for t in range(0, max_days + 1):
                    holder += "," + str(new_traj[bc][t])
                out_traj.write(holder + "\n")
            out_traj.close()
        else:
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

    # Writing outputs
    if doc:
        return barcodes, s, stats

    save_fitness(name + "_fitness.csv", barcodes, s, replicates, stats)
