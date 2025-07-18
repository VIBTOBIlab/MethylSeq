#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import warnings
from scipy.optimize import OptimizeWarning
import sys

def asymptotic_growth(x, beta0, beta1):
    return beta0 * np.arctan( beta1 * x ) #tanh cause it has an asymptotic growth similar to seq saturation


def derivative_asymptotic_growth(x, beta0, beta1):
    return beta0 * beta1 / (1 + (beta1 * x) ** 2)


def find_asymptote(params):
    beta0, beta1 = params
    asymptote = beta0 * np.pi / 2 # Limit for x that tends to Inf
    return asymptote

def plot_error_data(x_data,y_data,reads,output_plot,title,error):
    # Predict new data
    xpred = np.append(x_data, np.array( [1.2,1.4,1.6,1.8,2] ) )
    xpred_reads = np.array( reads*xpred )
    
    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the main data (on bottom x-axis)
    ax1.plot(x_data, y_data, marker="o", color="blue", linestyle="")

    # Plot the fit line
    ax1.set_xticks(xpred)

    # Set title and labels for the bottom x-axis
    ax1.set_title(title + "\n" + error, fontsize=12, color = "red", loc= "left")
    ax1.set_xlabel("Percentage of downsampling", fontsize=14)
    ax1.set_ylabel("Number of CpGs", fontsize=14)

    # Plot the second x-axis for xpred_reads (top x-axis)
    ax2 = ax1.twiny()  # Create a second x-axis sharing the same y-axis
    ax2.set_xlim(ax1.get_xlim())  
    ax2.set_xticks(xpred)
    ax2.tick_params(axis='x',pad=7)
    ax2.set_xticklabels([f'{x:.1e}' for x in xpred_reads], fontsize=10, rotation=45)
    ax2.set_xlabel("Number of reads", fontsize=14)
    # Customize the x-tick positions and move them to the top
    ax2.xaxis.set_ticks_position('top')
    ax1.xaxis.set_ticks_position('bottom')

    # Display grid and layout
    ax1.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()

    # Save the plot
    plt.savefig(output_plot)
    plt.close()
    print(f"Plot saved to {output_plot}")

def plot_data(x_data,y_data,reads,asymptote,params,output_plot,title):

    # Predict new data
    xpred = np.append(x_data, np.array( [1.2,1.4,1.6,1.8,2] ) )
    xpred_reads = np.array( reads*xpred )
    ypred = asymptotic_growth(xpred, *params)
    y_diff = [0,1]+[(ypred[i]/asymptote)*100 for i in range(2,len(ypred))]
    
    # Create the figure and axis
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the main data (on bottom x-axis)
    ax1.plot(x_data, y_data, marker="o", color="blue", linestyle="")

    # Plot the fit line
    ax1.plot(xpred, ypred, 'g-', label='fit: beta0=%5.3f, beta1=%5.3f' % tuple(params))
    ax1.plot(xpred, ypred, marker='|', color="black",linestyle="",markersize=10)
    ax1.set_xticks(xpred)

    # Add text labels for y_diff below the plot
    for i in range(4, len(xpred)):
        ax1.text(xpred[i], ypred[i] - 0.09 * ypred[i], f"{y_diff[i]:.1f}%", color="black", fontsize=10, ha="center")

    # Add asymptote line
    ax1.axhline(y=asymptote, color='grey', linestyle='--')

    # Set title and labels for the bottom x-axis
    ax1.set_title(title, fontsize=16, loc= "left")
    ax1.set_xlabel("Percentage of downsampling", fontsize=14)
    ax1.set_ylabel("Number of CpGs", fontsize=14)

    # Plot the second x-axis for xpred_reads (top x-axis)
    ax2 = ax1.twiny()  # Create a second x-axis sharing the same y-axis
    ax2.set_xlim(ax1.get_xlim())  
    ax2.set_xticks(xpred)
    ax2.tick_params(axis='x',pad=7)
    ax2.set_xticklabels([f'{x:.1e}' for x in xpred_reads], fontsize=10, rotation=45)
    ax2.set_xlabel("Number of reads", fontsize=14)
    # Customize the x-tick positions and move them to the top
    ax2.xaxis.set_ticks_position('top')
    ax1.xaxis.set_ticks_position('bottom')
    leg_patch = mpatches.Patch(label = r"   % : Sequencing saturation ($\frac{\hat y}{\text{asymptote}} \times 100$)")
    ax1.text(0, asymptote-asymptote*0.05, f"Asymptote = {asymptote:.2e}", color="grey", fontsize=12)
    plt.legend(handles=[leg_patch],loc="lower right",handletextpad=-1.0, handlelength=0)

    # Display grid and layout
    ax1.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()

    # Save the plot
    plt.savefig(output_plot)
    plt.close()
    print(f"Plot saved to {output_plot}")


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
    fit_success = True
    
    # If doesn't find optimal params, plot just the points without the curve_fit
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            # Fit the function to the data
            params, _ = curve_fit(asymptotic_growth, x_data, y_data, p0=[1, 1])
        asymptote = find_asymptote(params)
    except (RuntimeError, OptimizeWarning) as e:
        fit_success = False
        params = None
        asymptote = None
        error_msg = str(e)
        print(f"Curve fitting error: {e}")

    # Plot the data and save the plots
    title = list(set(data["sample"]))[0]
    reads = int(data["reads_counts"].tolist()[-1])

    if fit_success==True:
        plot_data(x_data,y_data,reads,asymptote,params,output_plot,title)
    else: plot_error_data(x_data,y_data,reads,output_plot,title,error_msg)


def select_sample(cpgs, reads,percentages):
    # Read the input file into a DataFrame
    column_names_cpgs = ["sample","percentage","min_counts","cpgs_counts"]
    column_names_reads = ["sample","percentage","reads_counts"]

    cpg_file = pd.read_csv(cpgs, sep=",",header=None,names=column_names_cpgs)
    reads_file = pd.read_csv(reads, sep=",",header=None,names=column_names_reads)
    data = pd.merge(cpg_file,reads_file,on=["sample","percentage"])
    samples = data["sample"].unique()
    for sample in samples:
        sample_data = data[data["sample"] == sample]
        min_counts = sample_data["min_counts"].unique()
        for min in min_counts:
            # Sort the dataframe to keep the ascending order for curve_fit function
            fin_df = sample_data[sample_data["min_counts"]==min].sort_values(by=['percentage'])
            output_plot = sample+"_"+str(min)+"x_plot.png"

            # If duplicated rows in the df, it means some processes were duplicated
            if fin_df.duplicated().any():
                print(f"Error: duplicated rows in {sample}: this might be due to error in Nextflow run (usually happening with -c vsc_ugent.config). Try to rerun without -resume flag.")
                sys.exit()

            plot_reads_vs_cpgs(fin_df,output_plot,percentages)


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
