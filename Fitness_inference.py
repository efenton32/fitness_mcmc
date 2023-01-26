# Imports
import sys
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

# Input parameters
population = sys.argv[1]
environment = sys.argv[2]
replicates = sys.argv[3]


