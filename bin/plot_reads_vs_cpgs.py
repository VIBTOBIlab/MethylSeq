import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from scipy.optimize import curve_fit
from scipy.optimize import fsolve


def asymptotic_growth(x, beta0, beta1):
    return beta0 * np.arctan( beta1 * x ) #tanh cause it has an asymptotic growth similar to seq saturation


def derivative_asymptotic_growth(x, beta0, beta1):
    return beta0 * beta1 / (1 + (beta1 * x) ** 2)


def plot_data(x_data,y_data,asymptote,params,output_plot,title):
    # Predict new data
    xpred = np.append(x_data, np.array( [1.2,1.4,1.6,1.8,2] ) )
    ypred = asymptotic_growth(xpred, *params)
    y_diff = [0,1]+[(ypred[i]-ypred[i-1])/ypred[i-1] for i in range(2,len(ypred))]
    
    # Plot the data
    plt.figure(figsize=(10, 6))
    plt.plot(x_data, y_data, marker="o", color="blue",linestyle="")
    plt.plot(xpred, ypred, 'g--',
         label='fit: beta0=%5.3f, beta1=%5.3f' % tuple(params))
    for i in range(4,len(xpred)):
        plt.text(xpred[i], ypred[i] - 0.07 * ypred[i], f"{y_diff[i]:.4f}", color="blue", fontsize=10, ha="center")
    plt.axhline(y=asymptote,color='r',linestyle='--')
    plt.text(0, asymptote - asymptote*0.1, f"Asymptote = {asymptote:.2e}", color="red", fontsize=12)
    plt.title(title, fontsize=16)
    plt.xticks(xpred, labels=[f"{x:.1f}" for x in xpred], fontsize=10)
    plt.xlabel("Percentage of downsampling", fontsize=14)
    plt.ylabel("Number of CpGs", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_plot)
    plt.close()
    print(f"Plot saved to {output_plot}")


def find_asymptote(params):
    beta0, beta1 = params
    asymptote = beta0 * np.pi / 2 # Limit for x that tends to Inf
    return asymptote


def plot_reads_vs_cpgs(data, output_plot,percentages):
    """
    Plots the percentage of subsample vs. the number of CpGs from the input file.

    Parameters:
        data (str): Dataframe containing reads and CpGs.
        output_plot (str): Path to save the output plot.
    """

    # Add 0 to the x and y data
    x_data = np.array( [0] + percentages )
    y_data = np.array( [0] + data["cpgs_counts"].tolist() )

    # Fit the function to the data
    params, _ = curve_fit(asymptotic_growth, x_data, y_data, p0=[1, 1])

    # Find the asymptote value
    asymptote = find_asymptote(params)
    
    # Plot the data and save the plots
    title = list(set(data["sample"]))[0]
    plot_data(x_data,y_data,asymptote,params,output_plot,title)


def select_sample(cpgs, reads,percentages):
    # Read the input file into a DataFrame
    column_names_cpgs = ["sample","percentage","cpgs_counts"]
    column_names_reads = ["sample","percentage","reads_counts"]

    cpg_file = pd.read_csv(cpgs, sep=",",header=None,names=column_names_cpgs)
    reads_file = pd.read_csv(reads, sep=",",header=None,names=column_names_reads)
    data = pd.merge(cpg_file,reads_file,on=["sample","percentage"])
    samples = data["sample"].unique()
    for sample in samples:
        sample_data = data[data["sample"] == sample]
        output_plot = sample+"_plot.png"
        plot_reads_vs_cpgs(sample_data,output_plot,percentages)


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Plot Reads vs. CpGs from two CSV files.")
    parser.add_argument("--cpgs_file", type=str, help="Path to the input CSV file.")
    parser.add_argument("--read_file", type=str, help="Path to the input CSV file.")
    parser.add_argument("--percentages",type=str,help="List of downsampling percentages.")

    args = parser.parse_args()
    percentages = [float(p) for p in args.percentages.split(",")]

    # Call the plotting function
    select_sample(args.cpgs_file,args.read_file,percentages)
