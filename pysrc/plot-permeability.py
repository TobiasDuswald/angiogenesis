import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

# Use latex in the figures
plt.rc("text", usetex=True)

# Get the directory of this script
import os

mydir = os.path.dirname(os.path.realpath(__file__))

# Define the filename
filename = os.path.join(mydir, "..", "data", "permeability_values.txt")

# Load the data
data = np.loadtxt(filename, delimiter=",")

# Compute the time axis, assuming that the data is sampled every 5 minutes
time = np.arange(0, len(data)) * 5 / 60 / 24 - 100

# Plot the data
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
plt.plot(time, 1 + data)
ax.set_xlabel("Time (Days)")
ax.set_ylabel(r"$\varphi (t) = 1 + \chi (t)$")

ax.axvspan(2, 3, alpha=0.2, color="fuchsia", label="TRA")
ax.axvspan(5, 6, alpha=0.2, color="fuchsia")
ax.axvspan(8, 9, alpha=0.2, color="lightseagreen", label="DOX")

# Define the ticks: 0, 2, 4, ..., 20
ticks = np.arange(0, 21, 2)

# Set the ticks
plt.xticks(ticks)

plt.xlim(0, 20)
plt.tight_layout()
plt.legend()
plt.show()
