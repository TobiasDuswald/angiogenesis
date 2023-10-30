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

import numpy as np
import matplotlib.pyplot as plt
import shutil, os

# Helper function h, smooth heaviside function
def h(x, a, b, xx, dt):
    return 1 - np.exp(-(a + 1 / (1 + np.exp(2 * b * (x - xx)))) * dt)


# Helper function l, linear increase
def l(x, c, xxx, dt):
    tmp = (x - xxx) / (1 - xxx)
    # Element wise maximum
    tmp = np.maximum(tmp, 0)
    return 1 - np.exp(-c * tmp * dt)


def plot(fn, func_name, x, dt, xx, a, b):
    # Define function
    if func_name == "h":
        f = h
    elif func_name == "l":
        f = l
    else:
        raise ValueError("Unknown function name")

    # Plot
    if func_name == "h":
        for aa in a:
            for bb in b:
                for xxx in xx:
                    plt.plot(
                        x,
                        f(x, aa, bb, xxx, dt),
                        label=r"$a={}, b={}, ".format(aa, bb)
                        + r"\bar{x}"
                        + r"={}$".format(xxx),
                    )
        plt.ylabel(r"$h(x;a,b,\bar{x})$")
    elif func_name == "l":
        for cc in a:
            for xxx in xx:
                plt.plot(
                    x,
                    f(x, cc, xxx, dt),
                    label=r"$c={}, ".format(cc)
                    + r"\bar{x}"
                    + r"={}$".format(xxx),
                )
        plt.ylabel(r"$l(x;c,\bar{x})$")
    for xxx in xx:
        # Plot a vertical, dotted, light gray line at x=xxx
        plt.axvline(x=xxx, color="lightgray", linestyle="dotted")
    plt.gcf().set_size_inches(6, 4)
    plt.legend(loc="best")
    plt.xlabel(r"$x$")
    plt.savefig("{}.pdf".format(fn))
    plt.savefig("{}.png".format(fn))
    plt.close()
    plt.clf()


def main():
    # Fixed parameters
    x = np.linspace(0, 1, 1000)
    dt = 0.01

    # Activate latex annotations for matplotlib
    plt.rc("text", usetex=True)

    # Get director of this script
    dir = os.path.dirname(os.path.realpath(__file__))

    # Create empty output director
    outdir = os.path.join(dir, "..", "results/helper_functions")
    # Absolute path
    outdir = os.path.abspath(outdir)
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Define parameters
    a = [1, 2, 3, 4, 5]
    b = [30]
    xx = [0.5]

    # Plot helper function h
    plot(os.path.join(outdir, "h_a"), "h", x, dt, xx, a, b)

    # Define parameters
    a = [2]
    b = [1, 5, 15, 30, 50]
    xx = [0.5]

    # Plot helper function h
    plot(os.path.join(outdir, "h_b"), "h", x, dt, xx, a, b)

    # Define parameters
    a = [2]
    b = [30]
    xx = [0.15, 0.3, 0.5, 0.7, 0.85]

    # Plot helper function h
    plot(os.path.join(outdir, "h_xbar"), "h", x, dt, xx, a, b)

    # Define parameters
    c = [2]
    xx = [0.15, 0.3, 0.5, 0.7, 0.85]

    # Plot helper function l
    plot(os.path.join(outdir, "l_xbar"), "l", x, dt, xx, c, [])

    # Define parameters
    c = [1, 2, 3, 4, 5]
    xx = [0.2]

    # Plot helper function l
    plot(os.path.join(outdir, "l_a"), "l", x, dt, xx, c, [])


if __name__ == "__main__":
    main()
