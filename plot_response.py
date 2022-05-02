#!/bin/python

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Create figure with correct size and 4 plot regions
    f = plt.figure(figsize=(6, 4))
    plt.title("Cellular Response function")

    attack = 0.5
    x0 = 5.0
    x_low = 0
    x_high = 20
    x = np.linspace(x_low,x_high)
    y = 0.5*(1 + special.erf(attack*x-x0))

    plt.xticks(
        [x_low, (x_low+x_high)/2, x_high],
        ["0", "$\\rho_0$", "$\\rho_{free}$"]
    )

    plt.plot(x, y, label="Response function", color="k")

    plt.legend()
    plt.xlabel("Substrate $\\rho$")
    plt.tight_layout()
    plt.savefig("ResponseFunction.png")
    plt.clf()