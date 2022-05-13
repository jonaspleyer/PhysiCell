#!/usr/bin/env python3

import glob
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    files = glob.glob("output/output*_microenvironment0.mat")
    data = [loadmat(file)["multiscale_microenvironment"]
            .reshape((6, 40, 40))
            for file in files]

    n = 0
    for i in range(2):
        a = np.average(data[n][i+4], axis=1)
        print(a)
        o = np.std(data[n][i+4], axis=1)
        plt.errorbar(
            range(len(a)),
            a,
            yerr=o,
            label="Reactant " + str(i)
        )
    # plt.colorbar(orientation='vertical')
    plt.yscale("log")
    plt.legend()
    plt.savefig("results_density_profile.png")