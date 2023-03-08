# This script takes the data from the rat tumor (diam, length) and tests
# different probability distributions to fit the data. The best fit is
# determined by the Kolmogorov-Smirnov test.

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import pandas as pd
import seaborn as sns

sns.set_style("whitegrid")

# Define all distributions to test
dist_names = [
    "alpha",
    "anglit",
    "arcsine",
    "argus",
    "beta",
    "betaprime",
    "bradford",
    "burr",
    "burr12",
    "cauchy",
    "chi",
    "chi2",
    "cosine",
    "crystalball",
    "dgamma",
    "dweibull",
    "expon",
    "exponnorm",
    "exponweib",
    "exponpow",
    "f",
    "fatiguelife",
    "fisk",
    "foldcauchy",
    "foldnorm",
    "genlogistic",
    "gennorm",
    "genpareto",
    "genexpon",
    "genextreme",
    "gausshyper",
    "gamma",
    "gengamma",
    "genhalflogistic",
    "genhyperbolic",
    "geninvgauss",
    "gompertz",
    "gumbel_r",
    "gumbel_l",
    "halfcauchy",
    "halflogistic",
    "halfnorm",
    "halfgennorm",
    "hypsecant",
    "invgamma",
    "invgauss",
    "invweibull",
    "johnsonsb",
    "johnsonsu",
    "kappa4",
    "kappa3",
    # "kstwo",
    "kstwobign",
    "laplace",
    "laplace_asymmetric",
    "levy",
    "levy_l",
    "logistic",
    "loggamma",
    "loglaplace",
    "lognorm",
    "loguniform",
    "lomax",
    "maxwell",
    "mielke",
    "moyal",
    "nakagami",
    "ncx2",
    "ncf",
    "nct",
    "norm",
    "norminvgauss",
    "pareto",
    "pearson3",
    "powerlaw",
    "powerlognorm",
    "powernorm",
    "rdist",
    "rayleigh",
    "rice",
    "recipinvgauss",
    "semicircular",
    "skewcauchy",
    "skewnorm",
    "t",
    "trapezoid",
    "triang",
    "truncexpon",
    "truncnorm",
    "tukeylambda",
    "uniform",
    "vonmises",
    "vonmises_line",
    "wald",
    "weibull_min",
    "weibull_max",
    "wrapcauchy",
]


def get_data():
    df = pd.read_csv("data/vessel-diameter.txt")
    # Strip the column names
    df.columns = df.columns.str.strip()
    return df

def compute_volume(df):
    # Get numpy arrays for diameter and length of cylinders
    diameters = df["diam"].to_numpy()
    lengths = df["length"].to_numpy()
    # Compute the volume of each cylinder
    volumes = np.pi * (diameters / 2) ** 2 * lengths
    # Compute the volume fraction
    volume_fraction = np.sum(volumes) / (550*550*230)
    # Print volume and volume fraction
    print("Volume: " + str(np.sum(volumes)))
    print("Volume fraction: " + str(volume_fraction))


def get_best_distribution(data):
    dist_results = []
    df_data = {"dist": [], "p": [], "params": []}
    params = {}
    for dist_name in dist_names:
        dist = getattr(st, dist_name)
        param = dist.fit(data)
        params[dist_name] = param
        # Applying the Kolmogorov-Smirnov test
        D, p = st.kstest(data, dist_name, args=param)
        print("p value for " + dist_name + " = " + str(p))
        dist_results.append((dist_name, p))
        df_data["dist"].append(dist_name)
        df_data["p"].append(p)
        df_data["params"].append(param)
    # select the best fitted distribution
    best_dist, best_p = max(dist_results, key=lambda item: item[1])
    # # store the name of the best fit and its p value
    print(
        "\033[92m" + "Best fitting distribution: " + str(best_dist) + "\033[0m"
    )
    print("\033[92m" + "Best p value: " + str(best_p) + "\033[0m")
    print(
        "\033[92m"
        + "Parameters for the best fit: "
        + str(params[best_dist])
        + "\033[0m"
    )
    # Create a dataframe with the results and sort it by p value
    df = pd.DataFrame(df_data)
    df = df.sort_values(by="p", ascending=False)
    print(df.head(10))
    print("...")
    print(df.tail(10))
    return best_dist, best_p, params[best_dist]

def plot_hist_and_fit(data, best_dist, best_p, params, type="diam"):
    # Plot for comparison
    plt.figure(figsize=(7, 5))
    plt.hist(data, bins="auto", density=True, alpha=0.5, label="Data")
    # Save the parameters used by the fit
    dist = getattr(st, best_dist)
    # Update the plot
    xmin = np.min(data)
    xmax = np.max(data)
    x = np.linspace(xmin, xmax, 1000)
    # Get PDF
    pdf = dist.pdf(x, *params)
    # Plot PDF
    plt.plot(x, pdf, label="PDF", linewidth=2)
    if type == "diam":
        plt.xlabel("Diameter (µm)")
    else:
        plt.xlabel("Length (µm)")
    plt.ylabel("Frequency")
    plt.legend(loc="best")
    # # Add small box box with the best fit and the p value
    # plt.text(
    #     0.95,
    #     0.95,
    #     "Best fit: " + str(best_dist) + "\np value: " + str(best_p)[0:6],
    #     fontsize=12,
    #     horizontalalignment="right",
    #     verticalalignment="top",
    #     transform=plt.gca().transAxes,
    # )
    # plt.show()
    # Save the plot as pdf
    plt.savefig(type + "-hist.pdf", bbox_inches="tight")



def main():
    data = get_data()

    # Print how many samples we have
    print("Number of samples: " + str(len(data)))

    # Compute the volume fraction
    compute_volume(data)

    # Describe the pandas data
    print(data.describe())

    # Compute the Pearson and Spearman correlation coefficients
    print("Pearson correlation coefficient: " + str(data.corr(method="pearson")))
    print("Spearman correlation coefficient: " + str(data.corr(method="spearman")))
    print("Kendall correlation coefficient: " + str(data.corr(method="kendall")))

    y = data["diam"]
    best_dist, best_p, params = get_best_distribution(y)
    plot_hist_and_fit(y, best_dist, best_p, params, "diam")

    y = data["length"]
    best_dist, best_p, params = get_best_distribution(y)
    plot_hist_and_fit(y, best_dist, best_p, params, "length")


if __name__ == "__main__":
    main()
