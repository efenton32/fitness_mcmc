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

out = "LTEE"
pop = "pop"

environment_list = ["aceH", "aceL", "ampL", "chlL", "cipL", "fruH", "fruL", "gluH", "gluL", "lacH", "lacL", "lbrH",
                    "lbrL", "malH", "malL", "rifL"]

for p in range(1, 5):
    for e in environment_list:
        print(p, e)
        #f.day_fitness(out, pop + str(p), e, point_fit=2)
        s, stats, barcodes = f.read_s(out, pop + str(p), e)
        f.fitness_plot(barcodes, s, stats, out + "_" + pop + str(p) + "_" + e + "_fit_day1.png")