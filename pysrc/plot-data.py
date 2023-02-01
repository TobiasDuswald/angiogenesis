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
        ax.plot(data.index, data["mean_{}".format(i)], label="mean")
        ax.fill_between(
            data.index,
            data["lower_{}".format(i)],
            data["upper_{}".format(i)],
            alpha=0.5,
            label="std",
        )

        # ax.set_ylim(-10, 510)

        # Set the labels
        ax.set_xlabel("Days")
        ax.set_ylabel("Tumor volume ($mm^3$)")
        # ax.set_title("Group {}".format(i))

        # Set the legend
        ax.legend()

        # Do not cut off the labels
        fig.tight_layout()

        # Save the figure
        fig.savefig(os.path.join(results_dir, "growth_group_{}.pdf".format(i)))

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
    results_dir = os.path.join(dir_path, "..", "results", "data-visualization")
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