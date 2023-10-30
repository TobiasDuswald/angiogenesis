# -----------------------------------------------------------------------------
#
# Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------

# This script solves a non-linear optimization problem to find the parameters a,
# b, gamma, xbar that minimize the error between data and the
# actual h function. The optimization problem is solved using the
# scipy.optimize.minimize function. The optimization problem is solved for
# different values of dt. The results are stored in the file
# pysrc/hfunction_optimizer_results.txt.
# The results are also plotted in the file pysrc/hfunction_optimizer_results.png.

from absl import app
from absl import flags
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize

FLAGS = flags.FLAGS
flags.DEFINE_float("xbar", 0.16, "xbar")
flags.DEFINE_float("dt", 1.0, "dt")

# Data to fit
# ymin = 0.5e-5
# ymax = 5e-5
# data = np.array(
#     [
#         [0.0, ymin],
#         [0.1, ymin],
#         [0.16, (ymin+ymax)/2], # xbar
#         [0.4, ymax],
#         [0.5, ymax],
#         [0.6, ymax],
#         [0.7, ymax],
#         [0.8, ymax],
#         [0.9, ymax],
#         [1.0, ymax],
#     ]
# )
ymin = 0
ymax = 1
data = np.array(
    [
        [0.0, ymin],
        [0.05, 0.5],
        [0.1, ymax],
        [0.16, ymax],  # xbar
        [0.4, ymax],
        [0.5, ymax],
        [0.6, ymax],
        [0.7, ymax],
        [0.8, ymax],
        [0.9, ymax],
        [1.0, ymax],
    ]
)


# Helper function h, smooth heaviside function
def h(x, a, b, gamma, xx, dt):
    return 1 - np.exp(-(a + 1 / (gamma + np.exp(2 * b * (x - xx)))) * dt)
    # return (1 - np.exp(-(a + 1 / (1.0 + np.exp(2 * b * (x - xx)))) * dt))/gamma
    # return (a + 1 / (gamma + np.exp(2 * b * (x - xx)))) * dt
    # return (a + b * np.tanh((x - xx) * gamma)) * dt
    # return 1- np.exp(-(a + b * np.tanh((x - xx) * gamma)) * dt)
    # return 1 - np.exp(-(a + 1 / (1.0 + np.exp(2 * b * (x - xx)))) * dt) * gamma
    # return a + b * (x - xx) + gamma * dt


def main(argv):
    # Get parameters
    xbar = FLAGS.xbar
    dt = FLAGS.dt

    # Define objective function
    def objective(x):
        """RMS between h(x) and data"""
        rms = 0.0
        for i in range(len(data)):
            rms += (h(data[i, 0], x[0], x[1], x[2], xbar, dt) - data[i, 1]) ** 2
        rms = np.sqrt(rms / len(data))
        return rms

    # Initial guess
    x0 = np.array([0.0, 10.0, 1.0])
    if data[0, 1] < data[-1, 1]:
        # correct function shape
        x0[1] = -1.0
    x_backup = x0.copy()

    # Solve optimization problem
    res = minimize(
        objective,
        x0,
        method="nelder-mead",
        options={"xtol": 1e-8, "disp": True},
    )

    # Print results
    print("xbar = {}".format(xbar))
    print("data = {}".format(data))
    print("x0 = {}".format(x_backup))
    print("x = {}".format(res.x))
    print("rms = {}".format(res.fun))

    # Plot results
    x = np.linspace(0, 1, 100)
    plt.plot(x, h(x, res.x[0], res.x[1], res.x[2], xbar, dt), label="h(x)")
    # plt.plot(
    #     x,
    #     h(x, x_backup[0], x_backup[1], x_backup[2], xbar, dt),
    #     label="h(x) initial guess",
    # )
    plt.plot(data[:, 0], data[:, 1], "o", label="data")
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("h(x)")
    plt.show()


if __name__ == "__main__":
    app.run(main)
