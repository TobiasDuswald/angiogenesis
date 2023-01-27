import pandas as pd
import numpy as np

"""
Specific script to combine the data from the file data/growth_data.csv.
usage: python pysrc/combine_data.py
"""


class Group:
    def __init__(self, n, mean, std):
        self.n = n
        self.mean = mean
        self.std = std


def read_data(filename, index_col=None):
    data = pd.read_csv(filename, index_col=index_col)
    return data


def combine_groups(group1, group2):
    """
    Combine two groups of data.
    """
    # Calculate the new mean
    new_mean = (group1.n * group1.mean + group2.n * group2.mean) / (
        group1.n + group2.n
    )
    # Calculate the new standard deviation (CoPilot suggestion)
    # new_std = np.sqrt(
    #     (
    #         group1.n * (group1.std**2 + group1.mean**2)
    #         + group2.n * (group2.std**2 + group2.mean**2)
    #     )
    #     / (group1.n + group2.n)
    #     - new_mean**2
    # )

    # StackExchange suggestion:
    # https://stats.stackexchange.com/questions/117741/adding-two-or-more-means-and-calculating-the-new-standard-deviation
    new_std = np.sqrt(
        (
            (group1.n - 1) * group1.std**2
            + (group2.n - 1) * group2.std**2
            + +(group1.n * group2.n)
            * (
                group1.mean**2
                + group2.mean**2
                - 2 * group1.mean * group2.mean
            )
            / (group1.n + group2.n)
        )
        / (group1.n + group2.n - 1)
    )

    # Calculate the new number of data points
    new_n = group1.n + group2.n

    return Group(new_n, new_mean, new_std)


if __name__ == "__main__":
    # Get directory of this file
    import os

    dir_path = os.path.dirname(os.path.realpath(__file__))

    # Define filename
    filename = os.path.join(dir_path, "..", "data", "growth_data.csv")
    filename = os.path.abspath(filename)

    # Define results directory (base path is the directory of filename)
    results_dir = os.path.dirname(filename)

    # Read data
    print("Reading data from {}...".format(filename))
    data = read_data(filename, index_col="days")

    new_data = {"days": [], "n": [], "mean_1": [], "std_1": []}

    # Iterate the rows of the data
    for i, row in data.iterrows():
        # Break the loop for that that contains treatment
        if i > 34:
            break

        # Create two groups of data
        group1 = Group(7, row["mean_1"], row["std_1"])
        group2 = Group(8, row["mean_2"], row["std_2"])
        group3 = Group(7, row["mean_3"], row["std_3"])
        group4 = Group(7, row["mean_4"], row["std_4"])
        group5 = Group(6, row["mean_5"], row["std_5"])
        group6 = Group(7, row["mean_6"], row["std_6"])

        # Combine the groups
        group = combine_groups(group1, group2)
        group = combine_groups(group, group3)
        group = combine_groups(group, group4)
        group = combine_groups(group, group5)
        group = combine_groups(group, group6)

        # Store the results
        new_data["days"].append(i)
        new_data["n"].append(group.n)
        new_data["mean_1"].append(group.mean)
        new_data["std_1"].append(group.std)

    # Create a new DataFrame
    new_data = pd.DataFrame(new_data)
    new_data = new_data.set_index("days")
    print(new_data)

    # Drop the column "n"
    new_data = new_data.drop("n", axis=1)

    # Save the data to a CSV file
    new_filename = os.path.join(results_dir, "growth_data_combined.csv")
    print("Saving data to {}...".format(new_filename))
    new_data.to_csv(new_filename)
