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

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
Visualize the data from the file data/growth_data.csv.
usage: python pysrc/plot-data.py

Note: can also plot the merged data, adjust data source in main and possibly
change the figure size.
"""


def read_data(filename, index_col=None):
    data = pd.read_csv(filename, index_col=index_col)
    return data


def trasform_data(data):
    """
    Determine the number of columns 2N. Then, from 1 to N, create the columns
    lower_1, lower_2, ..., lower_N, upper_1, upper_2, ..., upper_N based on the
    columns mean_1, mean_2, ..., mean_N. The lower_i and upper_i are computed
    adding and subtracting the standard deviation std_i to mean_i, respectively.
    Return N and the transformed data.
    """
    # Determine the number of columns N
    N = int(len(data.columns) / 2)

    # Create the columns lower_1, ..., lower_N; and upper_1, ..., upper_N
    for i in range(1, N + 1):
        data["lower_{}".format(i)] = (
            data["mean_{}".format(i)] - data["std_{}".format(i)]
        )
        data["upper_{}".format(i)] = (
            data["mean_{}".format(i)] + data["std_{}".format(i)]
        )

    return N, data


def create_folder(results_dir):
    """
    Create the folder results_dir if it does not exist.
    """
    import os

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        print("Folder {} already exists".format(results_dir))
        overwrite = input("Overwrite? (y/n) ")
        if overwrite == "y":
            pass
        else:
            print("Folder {} not overwritten".format(results_dir))
            import sys

            sys.exit()


def plot_data(N, data, results_dir):
    """
    For each i in 1, 2, ..., N, plot the column mean. Use the columns lower_i
    and upper_i to fill the area between the mean and the standard deviation.
    Save each plot separately in the folder results_dir.
    """
    # Set the style of the plot
    sns.set_style("whitegrid")

    # Create the folder results_dir if it does not exist
    create_folder(results_dir)

    # Hard code treatment times for vertical lines
    treatment = {
        "1": {"DOX": [], "TRA": [], "SAL": [35, 38, 39]},
        "2": {"DOX": [39], "TRA": [], "SAL": [35, 38]},
        "3": {"DOX": [], "TRA": [35, 38], "SAL": [39]},
        "4": {"DOX": [35], "TRA": [36, 39], "SAL": []},
        "5": {"DOX": [39], "TRA": [35, 38], "SAL": []},
        "6": {"DOX": [35, 38], "TRA": [35, 38], "SAL": []},
    }

    # Plot the data
    for i in range(1, N + 1):
        print("Plotting group {}...".format(i))
        # Create the figure
        fig, ax = plt.subplots()

        # Set the size of the figure
        fig.set_size_inches(4, 3)

        # Set resolution of the figure
        fig.set_dpi(500)

        # Plot the data
        ax.plot(
            data.index,
            data["mean_{}".format(i)],
            label="mean",
            marker="o",
            markersize=3,
            color="#307EC9",
        )
        ax.fill_between(
            data.index,
            data["lower_{}".format(i)],
            data["upper_{}".format(i)],
            color="#307EC9",
            alpha=0.2,
            label="std",
        )

        # ax.set_ylim(-10, 510)

        # Plot vertical lines for treatment times: DOX (red, dashed),
        # TRA (green, dotted), SAL (orange, dotted-dashed)
        dox_legend_cntr = 0
        tra_legend_cntr = 0
        sal_legend_cntr = 0
        for t in treatment[str(i)]["TRA"]:
            if tra_legend_cntr == 0:
                ax.axvline(t, color="fuchsia", linestyle="--", label="TRA")
                tra_legend_cntr += 1
            else:
                ax.axvline(t, color="fuchsia", linestyle="--")
        for t in treatment[str(i)]["DOX"]:
            if dox_legend_cntr == 0:
                ax.axvline(
                    t, color="lightseagreen", linestyle=":", label="DOX"
                )
                dox_legend_cntr += 1
            else:
                ax.axvline(t, color="lightseagreen", linestyle=":")
        for t in treatment[str(i)]["SAL"]:
            if sal_legend_cntr == 0:
                ax.axvline(t, color="orange", linestyle="-.", label="SAL")
                sal_legend_cntr += 1
            else:
                ax.axvline(t, color="orange", linestyle="-.")

        # Set the labels
        ax.set_xlabel("Time (Days)")
        ax.set_ylabel("Tumor volume ($mm^3$)")
        # ax.set_title("Group {}".format(i))

        # Set the legend
        ax.legend()

        # Set location of the legend to the upper left corner
        ax.legend(loc="upper left")
        # Set ylim to 0, 3000
        ax.set_ylim(0, 3000)

        # Do not cut off the labels
        fig.tight_layout()

        # Limit the y-axis to 0 and 3000
        ax.set_ylim(0, 3000)

        # Use transparent background
        fig.patch.set_alpha(0)

        # Save the figure
        fig.savefig(os.path.join(results_dir, "growth_group_{}.png".format(i)))

        # Close the figure
        plt.close(fig)


if __name__ == "__main__":
    # Get directory of this file
    import os

    dir_path = os.path.dirname(os.path.realpath(__file__))

    # Define filename
    filename = os.path.join(dir_path, "..", "data", "growth_data.csv")
    # filename = os.path.join(dir_path, "..", "data", "growth_data_combined.csv")
    filename = os.path.abspath(filename)

    # Define results directory
    results_dir = os.path.join(
        dir_path, "..", "results", "data-visualization-2"
    )
    results_dir = os.path.abspath(results_dir)

    # Read data
    print("Reading data from {}...".format(filename))
    data = read_data(filename, index_col="days")

    # Transform data
    print("Transforming data...")
    N, data = trasform_data(data)

    # Plot data
    print("Plotting data...")
    plot_data(N, data, results_dir)
