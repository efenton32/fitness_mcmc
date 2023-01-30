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
pop = "pop1"

environment_list = ["aceH", "aceL", "ampL", "chlL", "cipL", "fruH", "fruL", "gluH", "gluL", "lacH", "lacL", "lbrH",
                    "lbrL", "malH", "malL", "rifL"]

for e in environment_list:
    print(e)
    f.fitness_pipeline(out, pop, e, format=True)