import numpy as np
from absl import app
from absl import flags

FLAGS = flags.FLAGS

flags.DEFINE_float("r", "1", "The radius of the agent.")
flags.DEFINE_float("ra", "1.1", "The action radius of an agent.")
flags.DEFINE_float("rn", ".5", "The nuclear radius of an agent.")
flags.DEFINE_float("ca", "0.1", "Coefficient for adhesive force]")
flags.DEFINE_float(
    "cr", "1", "Coefficient for repulsive force"
)

def AdhesiveForce(d):
  r = 2 * FLAGS.r
  ra = 2 * FLAGS.ra
  rn = 2 * FLAGS.rn
  if (0 < d < ra):
    factor = (d/ra - 1)**2
  else:
    factor = 0
  return -factor

def RepulsiveForce(d):
  r = 2 * FLAGS.r
  ra = 2 * FLAGS.ra
  rn = 2 * FLAGS.rn
  if (0 < d < rn):
    factor = -(rn*d/(r*r)-2*d/r+1)
  elif (rn <= d <= r ):
    factor = -(d*d/(r*r)-2*d/r+1)
  else:
    factor = 0
  return - factor


def main(argv):
  max_distance = FLAGS.ra * 2.1
  distances = np.arange(0.001, max_distance, 0.01)
  adhesive_forces = np.array([AdhesiveForce(d) for d in distances])
  repulsive_forces = np.array([RepulsiveForce(d) for d in distances])
  total_forces = FLAGS.ca * adhesive_forces + FLAGS.cr * repulsive_forces

  import matplotlib.pyplot as plt
  # Use latex annotations
  plt.rc('text', usetex=True)
  plt.plot(distances, adhesive_forces, label=r"$F_a$")
  plt.plot(distances, repulsive_forces, label=r"$F_r$")
  plt.plot(distances, total_forces, label=r"$c_a \cdot F_a + c_r \cdot F_r$")
  # plot a vertical line at the action radius
  plt.axvline(x=2*FLAGS.ra, color="black", linestyle="-.", label=r"$2 \cdot r_a$")
  plt.axvline(x=2*FLAGS.r, color="black", linestyle="--", label=r"$2 \cdot r$")
  plt.axvline(x=2*FLAGS.rn, color="black", linestyle="dotted", label=r"$2 \cdot r_n$")
  plt.xlabel("Distance")
  plt.ylabel("Force")
  plt.legend()
  plt.grid()
  plt.show()


if __name__ == "__main__":
    app.run(main)
