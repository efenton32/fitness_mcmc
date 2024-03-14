# Standard Imports
import numpy as np
import pymc3 as pm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import arviz as az

# fitness_mcmc package imports
import fitness_inference as f
import fitness_mcmc
import fitness_mcmc.data_io as io
import fitness_mcmc.fitness_mcmc as m
import fitness_mcmc.format_raw as raw

out = "MEGA"

#env = ["37C", "YPD", "YPA"]
#gen = ["50", "150", "300", "550", "800", "1000"]
env = ["37C", "YPD", "YPA"]
gen = ["300", "1000"]

#gen = ["300", "550", "800", "1000"]
testing = False

if testing:
    f.fitness_pipeline(output="UMI", population="Gen" + gen[0], environment=env[0], assay=1, replicates=1, max_days=3, gens_per_day=10)
else:
    for g in gen:
        for e in env:
            for a in range(1, 19):
                print(g, e, a)
                if g == "300" and e == "YPA":
                    pass
                elif g == "1000" and e == "37C":
                    pass
                elif g == "1000" and e == "YPD" and a < 7:
                    pass
                else:
                    f.fitness_pipeline(output="UMI", population="Gen" + g, environment=e, assay=a, replicates=1, max_days=3, gens_per_day=10)

                #s, stats, barcodes = f.read_s(out, "Gen" + g, e)
                #f.fitness_plot(barcodes, s, stats, out + "_" + "Gen" + g + "_" + e + "_assay" + str(a) + "_fitComp.png")
