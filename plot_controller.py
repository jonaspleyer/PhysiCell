#!/bin/python

import numpy as np
import matplotlib.pyplot as plt

def get_acceptable_data(values, buffer=5, sigma=0.0001):
    acc = []
    for i, val in enumerate(values):
        subrange = values[max(0,i-buffer): min(len(values), i+buffer)]
        acc.append(abs(np.median(subrange)-val) <= np.std(subrange) / sigma)
    return acc

if __name__ == "__main__":
    for l in range(2):
        filename = "controller_logs_" + str(l+1) + ".txt"
        names = np.loadtxt(filename, dtype=str, delimiter=",", max_rows=1)
        results = np.loadtxt(filename, delimiter=",", skiprows=1)

        # Create figure with correct size and 4 plot regions
        f = plt.figure(figsize=(12, 8))
        f, axes = plt.subplots(4, 1, sharex=True)
        axes[0].set_title("Controller Metrics")
        for i, name in enumerate(names[1:-2]):
            acc = get_acceptable_data(results[:,i+1])
            axes[i].plot(results[:,0][acc], results[:,i+1][acc], label=name, color="k")
            axes[i].legend()
            # axes[-1].set_ylabel("Controller Module Responses")
            axes[i].tick_params(axis="x",direction="in")
        
        acc = get_acceptable_data(results[:,-1])
        axes[-1].plot(results[:,0][acc], results[:,-1][acc], label="Difference $\Delta=V_t-V_o$", color="k")
        axes[-1].legend()
        axes[-1].set_xlabel("Update steps [10 min]")
        f.tight_layout()
        f.subplots_adjust(hspace=0)
        f.savefig("results_" + str(l+1) + ".png")
        plt.clf()