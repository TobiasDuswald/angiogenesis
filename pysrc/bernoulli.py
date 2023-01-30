# This scipts plots the probability of a Bernoulli process with a given
# probability p and a given number of steps.
# usage:
# python pysrc/bernoulli.py --pmin -6 --pmax -2.5


from absl import app
from absl import flags
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

FLAGS = flags.FLAGS
flags.DEFINE_float("pmin", -5, "log10(pmin) (e.g. -5 for pmin=10^{-5}")
flags.DEFINE_float("pmax", -2.5, "log10(pmax) (e.g. -3 for pmax=10^{-3}")


def bernoulli(x, steps):
    return (1 - x) ** steps


def main(argv):
    # Print parameters
    pmin = FLAGS.pmin
    pmax = FLAGS.pmax
    print("pmin = 10^{}".format(pmin))
    print("pmax = 10^{}".format(pmax))

    # Define steps in minutes
    hour = 60
    quaterday = 6 * hour
    halfday = 12 * hour
    day = 2 * halfday

    # Plot
    x = np.logspace(pmin, pmax, 1000)
    sns.set()
    plt.rc("text", usetex=True)
    plt.plot(x, bernoulli(x, hour), label=r"1h ($N=60$)")
    plt.plot(x, bernoulli(x, quaterday), label=r"6h ($N=360$)")
    plt.plot(x, bernoulli(x, halfday), label=r"12h ($N=720$)")
    plt.plot(x, bernoulli(x, day), label=r"24h ($N=1440$)")
    plt.legend()
    plt.xlabel(r"$p$")
    plt.ylabel(r"$(1-p)^N$")
    plt.show()


if __name__ == "__main__":
    app.run(main)
