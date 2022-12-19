import pandas as pd
import numpy as np
import os
import sys
sys.path.append('..')
from fitness_mcmc import create_trajectories, sample_lineages

def _get_file_path(filename, data_dir):
    """
    Takes the file name and returns the absolute data file path given that
    the data file is found in data_dir

    Params:
        file_name[str]: file name of data file
        data_dir[str]: the directory where the file is located

    Returns:
        data_path[str]: the absolute path to file_name
    """
    start = os.path.abspath(__file__)
    start_dir = os.path.dirname(start)
    data_dir = os.path.join(start_dir, data_dir)
    data_path = os.path.join(start_dir, data_dir, filename)

    return data_path

def load_data(filename, data_dir = "experimental_data", load_metadata = False,
    return_ordered = False, delimiter = "\t"):
    """
    Takes data file path and outputs time generations and counts as numpy arrays

    Params:
        data_file[str]: Path to data_file
        load_metada [bool]: Load true s and f0 values from simulated data
        return_ordered [bool]: Automatically reorder lineages from most to least
            total counts

    Returns:
        data [pandas_dataframe]: Actual data where time and counts is stored
        time [array_like]: Time generations when sampling would occur
        counts [array_like]: Counts of genotypes
        frequencies [array_like]: Normalized counts of genotypes
    """
    data_file = _get_file_path(filename, data_dir)
    data = pd.read_csv(data_file, delimiter = delimiter)
    time = [float(i) for i in data.columns[1:]]
    counts = data.loc[:,data.columns[1:]].to_numpy()

    if load_metadata:
        metadata_file = _get_file_path(
            filename.split(".txt")[0] + "_metadata.txt", data_dir
        )
        metadata = pd.read_csv(metadata_file, delimiter = "\t").to_numpy()[:,1:]
        s_vals = metadata[:, 0]
        f0_vals = metadata[:, 1]

    if return_ordered:
        idx = np.flipud(np.argsort(np.sum(counts, axis = 1)))
        counts = counts[idx, :].astype("float")
        data = data.reindex(idx)

    if load_metadata:
        if return_ordered:
            s_vals = s_vals[idx]
            f0_vals = f0_vals[idx]
        return data, time, counts, s_vals, f0_vals
    else:
        return data, time, counts

def write_simulated_datafile(filename, N = 40, times = [5, 10, 25, 40, 45],
        s_range = 0.1, depth = 1000, s_vals = [], f0_vals = []):
    """
    Creates a textfile of simulated trajectories formated like a real datafile

    Params:
        filename [str]: Name of the output file.
        N [int]: Population size, i.e. number of genotypes. Automatically
            assumed if f0_vals or s_vals are included.
        times [array_like]: Times, in generations, to sample lineages.
        s_range [float]: Range of fitness values. Ignored if s_vals is included.
        depth [int_or_float]: Simulated read depth, affects noise.
        s_vals [array_like]: Fitness values for the population, optional.
        f0_vals [array_like]: Starting frequencies of the population, optional.
    """
    if len(f0_vals) > 0 or len(s_vals) > 0:
        if len(f0_vals) > 0 and len(s_vals) > 0 and len(f0_vals) != len(s_vals):
            raise ValueError("s_vals and f0_vals must have the same length.")
        N = max(len(f0_vals), len(s_vals))
    if len(f0_vals) == 0:
        f0_vals = np.random.random(N)
    if len(s_vals) == 0:
        s_vals = np.random.random(N) * s_range
    times = np.array(times)

    trajectory = create_trajectories(f0_vals, s_vals, times)
    sampled = pd.DataFrame(sample_lineages(trajectory, depth * N),
                           columns = times)
    metadata = pd.DataFrame({"s_vals": s_vals, "f0_vals": f0_vals})

    if ".txt" in filename:
        filename = filename.split(".txt")[0]

    sampled.to_csv(filename + ".txt", sep="\t", index_label = "BC")
    metadata.to_csv(filename + "_metadata.txt", sep="\t", index_label = "BC")
