# This script allows to plot the h function from the paper.
# usage:
# python pysrc/hfunction.py --xbar 0.1614 --a 0.000068 --b 50 --gamma 0.000408 --dt 1

from absl import app
from absl import flags
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

FLAGS = flags.FLAGS
flags.DEFINE_float("xbar", 0.5, "xbar")
flags.DEFINE_float("a", 0.1, "a")
flags.DEFINE_float("b", 0.1, "b")
flags.DEFINE_float("gamma", 1, "gamma")
flags.DEFINE_float("dt", 0.1, "dt")
flags.DEFINE_boolean("approx", False, "Plot approximation?")

# Helper function h, smooth heaviside function
def h(x, a, b, gamma, xx, dt):
    return 1 - np.exp(-(a + 1 / (gamma + np.exp(2 * b * (x - xx)))) * dt)


# Approximation of h
def h_approx(x, a, b, gamma, xx, dt):
    return (a + 1 / (gamma + np.exp(2 * b * (x - xx)))) * dt


def main(argv):
    # Print parameters
    dt = FLAGS.dt
    xxx = FLAGS.xbar
    aa = FLAGS.a
    bb = FLAGS.b
    gamma = FLAGS.gamma
    plot_approximation = FLAGS.approx

    print("xbar = {}".format(xxx))
    print("a = {}".format(aa))
    print("b = {}".format(bb))
    print("gamma = {}".format(gamma))
    print("dt = {}".format(dt))

    # Plot
    x = np.linspace(0, 1, 500)
    sns.set()
    plt.rc("text", usetex=True)
    plt.plot(x, h(x, aa, bb, gamma, xxx, dt), label=r"$h(x)$")
    if plot_approximation:
        plt.plot(
            x, h_approx(x, aa, bb, gamma, xxx, dt), label=r"$h_{approx}(x)$"
        )
        plt.legend()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$h(x;a,b,\bar{x})$")
    plt.show()


if __name__ == "__main__":
    app.run(main)
