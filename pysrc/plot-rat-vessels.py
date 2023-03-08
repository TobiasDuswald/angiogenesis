# This script plots the rat vessels in 3D. We have a list of start and end
# points, one for each vessel. We plot the connections between the points
# in 3D.

import numpy as np
import matplotlib.pyplot as plt

def parse_data(filename, cols_to_drop):
    with open(filename, 'r') as f:
        data = f.readlines()
        # Drop the first two lines
        data = data[2:]
        # Remove the first cols_to_drop characters
        data = [line[cols_to_drop:] for line in data]
        # Replace an arbitraty number of spaces with a single space
        import re
        data = [re.sub(' +', ' ', line) for line in data]
        # Strip all whitespaces from the beginning and end of the lines
        data = [line.strip() for line in data]
        # Split the lines into lists of strings
        data = [line.split(' ') for line in data]
        # Filter out the empty strings
        data = [[x for x in line if x != ''] for line in data]
        # Convert the strings to floats
        data = [[float(x) for x in line] for line in data]
        # Convert the list of lists to a numpy array
        data = [np.array(line) for line in data]
    start = []
    end = []
    for line in data:
        start.append(line[0:3])
        end.append(line[3:])
    # Print the first two and last two points
    print("Line 1: start = {}, end = {}".format(start[0], end[0]))
    print("Line 2: start = {}, end = {}".format(start[1], end[1]))
    print("...")
    print("Line -2: start = {}, end = {}".format(start[-2], end[-2]))
    print("Line -1: start = {}, end = {}".format(start[-1], end[-1]))
    return start, end

# Add a straight line connecting two points to the plot
def add_line_to_plot(ax, start, end, color):
    ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=color)

def get_metadata(usecase):
    metadata = {}
    if usecase == "rattumor":
        metadata["filename"] = "data/rattum98_0.txt"
        metadata["Title"] = "Rat Tumor"
        metadata["cols_to_drop"] = 38
    elif usecase == "ratbrain":
        metadata["filename"] = "data/brain99.txt"
        metadata["Title"] = "Rat Brain"
        metadata["cols_to_drop"] = 33
    else:
        raise ValueError("Invalid usecase")
    return metadata

def plot_histograms(start, end):
    length = []
    for i in range(len(start)):
        length.append(np.linalg.norm(end[i] - start[i]))
    # Plot normalized histogram, automatically determining the bins, and plot
    # the KDE as well
    import seaborn as sns
    sns.distplot(length, hist=True, kde=True,
                    bins=int(180/5), color = 'darkblue',
                    hist_kws={'edgecolor':'black'},
                    kde_kws={'linewidth': 4})
    plt.title('Histogram of vessel lengths')
    plt.xlabel('Length')
    plt.ylabel('Density')
    plt.show()

    # Save the length data to a text file
    with open("data/vessel-lengths.txt", 'w') as f:
        for l in length:
            f.write("{}\n".format(l))
    


def main():
    usecase = "rattumor"

    # Get the metadata for the usecase
    metadata = get_metadata(usecase)

    # Get the start and end points for the lines
    # start, end = get_data()
    start, end = parse_data(metadata["filename"], metadata["cols_to_drop"])

    # Create a figure and 3D axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    n_points = len(start)

    # Add the lines to the plot
    for i in range(n_points):
        add_line_to_plot(ax, start[i], end[i], 'r')
    
    plt.title(metadata["Title"])
    plt.show()

    # Plot the histograms of the of the vessels
    plot_histograms(start, end)

if __name__ == '__main__':
    main()
