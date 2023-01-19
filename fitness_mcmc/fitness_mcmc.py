import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import arviz as az
import csv

class Fitness_Model:
    def __init__(self, data, times=None,
                 reps=1, s_ref=0, prior="gauss", s_prior=None):
        """
        Initializes the Fitness_Model class

        Params:
            data [array-like]: lineage counts over time. Shape:
                (# lineages, # times)
            times [array-like]: times, in generations, where lineages were
                sampled
            reps [int]: number of replicates experiments
            s_ref [float or int]: fitness value for the reference lineage;
                defaults to 0
            prior [str]: Prior to choose for fitness values. Options are
                "gauss" [default], "flat", or "values"
            s_prior [array_like]: mean and standard deviation for priors on
                fitness, eg. from another experiment. Should have shape
                (# lineages - 1, 2)
        """
        if times is None:
            self.times = np.array([7, 14, 28, 42, 49]).reshape([1, -1])
        else:
            self.times = np.array(times).reshape([1, -1, 1])
        self.num_times = len(self.times[0,:])
        self.reps = reps
        self.data = np.array(data).reshape([-1, self.num_times, reps])
        self.N = len(data[:, 0])
        self.model = pm.Model()
        self.s_ref_val = s_ref
        self.prior = prior
        if prior == "values":
            self.s_prior = s_prior

        with self.model:
            self.s_ref = pm.math.constant(s_ref, ndim = 3)

            if self.prior == "gauss":
                self.mu = pm.Uniform("mu", -0.5, 0.2)
                self.sigma = pm.Uniform("sigma", 0.01, 1.0)
                self.s = pm.Normal("s", self.mu, self.sigma,
                    shape = (self.N - 1, 1, 1)
                )
            elif self.priors == "flat":
                self.s = pm.Flat("s", shape = (self.N - 1, 1, 1))
            elif self.priors == "values":
                self.s = pm.Normal("s", self.s_prior[:,0], self.s_prior[:,1],
                    shape = (self.N - 1, 1, 1)
                )

            self.s_tot = pm.math.concatenate((self.s_ref, self.s))
            self.f0_list = [pm.Dirichlet(f"f0_{rep}", a = np.ones(self.N)).reshape((-1, 1, 1))
                            for rep in range(self.reps)]

            self.f0 = pm.math.concatenate(self.f0_list, axis = 2)

            self.f_tot = (self.f0 * pm.math.exp(self.s_tot * self.times)
                / pm.math.sum(self.f0 * pm.math.exp(self.s_tot * self.times),
                axis = 0)
            )

            self.n_obs = pm.Multinomial("n_obs",
                np.sum(self.data, axis = 0).reshape((1, self.num_times, self.reps)).transpose((1,2,0)),
                p = self.f_tot.transpose((1,2,0)), observed = self.data.transpose((1,2,0))
            )

    def mcmc_sample(self, draws, tune = 4000, **kwargs):
        """
        Markov-chain Monte Carlo sample the posterior distribution of the
        pymc3 model
        """
        with self.model:
            self.trace = pm.sample(draws,
                                   tune = tune,
                                   return_inferencedata = True,
                                   init = "adapt_diag",
                                   **kwargs)

    def plot_mcmc_posterior(self, save_path = None):
        """
        Plots the posterior for several MCMC sampled parameters
        """
        with self.model:
            if not save_path:
                az.plot_posterior(self.trace)
            else:
                axes = az.plot_posterior(self.trace)
                fig = axes.ravel()[0].figure
                fig.savefig(save_path + "posterior.png")

    def plot_mcmc_trace(self, save_path = None):
        """
        Plots the trace of a sampled MCMC posterior distribution
        """
        with self.model:
            if not save_path:
                az.plot_trace(self.trace)
            else:
                axes = az.plot_trace(self.trace)
                fig = axes.ravel()[0].figure
                fig.savefig(save_path + "trace.png")

    def find_MAP(self, **kwargs):
        """
        Finds the MAP estimate for lineage fitnesses and starting frequencies
        """
        self.map_estimate = pm.find_MAP(model = self.model, **kwargs)

    def save_MAP(self, save_file, barcodes=None, **kwargs):
        """
        Saves the MAP estimate to a csv file

        Parameters:
            save_file [str]: output file to save to
            barcodes [list]: barcodes or other unique identifier for lineages
        """
        if barcodes is None:
            barcodes = np.arange(len(self.map_estimate["s"])+1)
        with open(save_file, "w") as sf:
            writer = csv.writer(sf, delimiter="\t")
            writer.writerow([barcodes[0], self.s_ref_val])
            for bc, s_val in zip(barcodes[1:], self.map_estimate["s"]):
                writer.writerow([bc, s_val[0][0]])

    def plot_MAP_estimate(self, type="log_y"):
        """
        Plots lineage trajectories from the MAP estimate

        Parameters:
            type [str]: either "log_y" or "lin", sets the y axis scale
        """
        self.f_pred = np.zeros_like(self.data)
        self.f_pred[1:, :] = (self.map_estimate["f0"][1:, None]
            * np.exp(self.map_estimate["s"] * self.times)
        )
        self.f_pred[0] = (self.map_estimate["f0"][0]
            * np.exp(self.s_ref_val * self.times) )
        self.f_pred /= np.sum(self.f_pred, axis = 0)

        fig, axs = plt.subplots(1,2)
        if type == "log_y":
            axs[0].semilogy(self.times.T,
                            (self.data / np.sum(self.data, axis = 0)).T
            )
            axs[1].semilogy(self.times.T, self.f_pred.T)
        elif type == "lin":
            axs[0].plot(self.times.T,
                        (self.data / np.sum(self.data, axis = 0)).T
            )
            axs[1].plot(self.times.T, self.f_pred.T)
        axs[0].set_xlabel("Generations")
        axs[1].set_xlabel("Generations")
        axs[0].set_ylabel("Lineage frequency")
        axs[0].set_title("Data")
        axs[1].set_title("Reconstructed")
        plt.show()

def normalize_func(x):
    """
    Normalizes lineage frequencies to sum to 1

    Parameters:
        x [array_like]: frequency vector to be normalized
    """
    return x / np.sum(x, axis = 0)

def create_trajectories(f0, s, times, normalize = True):
    """
    Simulates lineage trajectores given initial frequency and fitnesses

    Parameters:
        f0 [array-like]: initial lineage frequencies
        s [array-like]: lineage fitnesses
        times [array-like]: times, in generations, to sample lineage frequencies
        normalize [bool]: if True, normalizes lineage frequencies at each
            sampling time
    Returns:
        f_traj [numpy array]: array of lineage frequencies sampled at times
            given by "times"
    """
    f0 = np.array(f0).reshape([len(f0), -1])
    s = np.array(s).reshape([len(s), -1])
    times = np.array(times).reshape([-1, len(times)])
    f_traj = f0 * np.exp(s * times)
    if normalize:
        f_traj = normalize_func(f_traj)
    return f_traj

def sample_lineages(f, num_samples):
    """
    Returns lineage counts Poisson sampled from their true frequencies

    Parameters:
        f [array_like]: true lineage frequencies
        num_samples [int]: total number of samples to draw, should be order
            100 * num_lineages
    Returns:
        n_sampled [numpy array]: number of samples measured from each lineage
    """
    n_expected = f * num_samples
    n_sampled = np.random.poisson(n_expected)
    return n_sampled
