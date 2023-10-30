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

# This scipts plots the effects of the choice of the timestep, e.g. the error
# of taking larger time steps.
# usage:
# python pysrc/timestep.py


from absl import app
from absl import flags
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

FLAGS = flags.FLAGS
flags.DEFINE_float("dtmin", -6, "log10(dtmin) (e.g. -5 for dtmin=10^{-5}")
flags.DEFINE_float("dtmax", 1, "log10(dtmax) (e.g. -3 for dtmax=10^{-3}")


def bernoulli(dt, rate, T):
    return (1 - rate * dt) ** (T / dt)


def main(argv):
    # Print parameters
    dtmin = FLAGS.dtmin
    dtmax = FLAGS.dtmax
    print("dtmin = 10^{}".format(dtmin))
    print("dtmax = 10^{}".format(dtmax))

    # Define steps in minutes
    hour = 60
    quaterday = 6 * hour
    halfday = 12 * hour
    day = 2 * halfday

    # Plot
    x = np.logspace(dtmin, dtmax, 1000)
    sns.set()
    plt.rc("text", usetex=True)

    rate = 0.0001
    plt.plot(
        x,
        bernoulli(x, rate, halfday) / bernoulli(x[0], rate, halfday),
        label=r"$T=12h$" + r", $r=10^{-4} min^{-1}$",
    )
    plt.plot(
        x,
        bernoulli(x, rate, day) / bernoulli(x[0], rate, day),
        label=r"$T=24h$" + r", $r=10^{-4}min^{-1}$",
    )
    rate = 0.0005
    plt.plot(
        x,
        bernoulli(x, rate, halfday) / bernoulli(x[0], rate, halfday),
        label=r"$T=12h$" + r", $r=5 \cdot 10^{-4}min^{-1}$",
    )
    plt.plot(
        x,
        bernoulli(x, rate, day) / bernoulli(x[0], rate, day),
        label=r"$T=24h$" + r", $r=5 \cdot 10^{-4}min^{-1}$",
    )
    rate = 0.001
    plt.plot(
        x,
        bernoulli(x, rate, halfday) / bernoulli(x[0], rate, halfday),
        label=r"$T=12h$" + r", $r=10^{-3}min^{-1}$",
    )
    plt.plot(
        x,
        bernoulli(x, rate, day) / bernoulli(x[0], rate, day),
        label=r"$T=24h$" + r", $r=10^{-3}min^{-1}$",
    )
    rate = 0.002
    plt.plot(
        x,
        bernoulli(x, rate, halfday) / bernoulli(x[0], rate, halfday),
        label=r"$T=12h$" + r", $r=2 \cdot 10^{-3}min^{-1}$",
    )
    plt.plot(
        x,
        bernoulli(x, rate, day) / bernoulli(x[0], rate, day),
        label=r"$T=24h$" + r", $r=2 \cdot 10^{-3}min^{-1}$",
    )
    plt.xscale("log")
    plt.legend()
    plt.xlabel(r"$\Delta t [min]$")
    plt.ylabel(
        r"$(1-r \cdot \Delta t)^{T/\Delta t} / (1-r \cdot \Delta t_{\min})^{T/\Delta t_{\min}} $"
    )
    plt.show()


if __name__ == "__main__":
    app.run(main)
