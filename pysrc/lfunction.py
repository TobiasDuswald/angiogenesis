# This scipts plots the l function from the paper.
# usage:
# python pysrc/lfunction.py --xbar 0.2 --c 3 --dt 0.1


from absl import app
from absl import flags
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

FLAGS = flags.FLAGS
flags.DEFINE_float("xbar", 0.5, "xbar")
flags.DEFINE_float("c", 0.1, "c")
flags.DEFINE_float("dt", 0.1, "dt")

# Helper function l, linear increase
def l(x, c, xxx, dt):
    tmp = (x - xxx) / (1 - xxx)
    # Element wise maximum
    tmp = np.maximum(tmp, 0)
    return 1 - np.exp(-c * tmp * dt)


def main(argv):
    # Print parameters
    dt = FLAGS.dt
    xxx = FLAGS.xbar
    cc = FLAGS.c
    print("xbar = {}".format(xxx))
    print("c = {}".format(cc))
    print("dt = {}".format(dt))

    # Plot
    x = np.linspace(0, 1, 500)
    sns.set()
    plt.rc("text", usetex=True)
    plt.plot(x, l(x, cc, xxx, dt))
    plt.xlabel(r"$x$")
    plt.ylabel(r"$l(x;c,\bar{x})$")
    plt.show()


if __name__ == "__main__":
    app.run(main)
