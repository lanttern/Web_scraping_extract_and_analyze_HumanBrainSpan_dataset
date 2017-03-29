from pylab import *
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

def heat_map(file1):
    """this function generate heatmap from dataframe"""
    df = pd.read_csv(file1, index_col=0)
    #print df
    df = df.div(df.mean(axis=1), axis=0)
    df = df.fillna(0)
    #print df
    #Plot it out
    fig, ax = plt.subplots()
    hmap = ax.pcolor(df, cmap = plt.cm.RdYlBu_r, alpha=1)
    plt.colorbar(hmap)
    #Format
    fig = plt.gcf()
    fig.set_size_inches(6,12)
    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(df.shape[0]) + 0.5, minor=False)
    #ax.set_xticks(np.arange(df.shape[1]) + 0.5, minor=False)
    
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    #fig.colorbar()
    # Set the labels
    #labels = ["ITC","STC","V1C","PFC", "CB"]
    #ax.set_xticklabels(labels, fontsize = 12)
    #ax.set_yticklabels(df.index, fontsize = 10)
    
    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

# Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
"""
    plt.yticks(np.arange(np.arange(df.shape[0]) + 0.5, len(df.index)), df.index)
    plt.xticks(np.arange(np.arange(df.shape[1]) + 0.5, len(df.column)), df.column)
    plt.tight_layout()
    plt.colorbar()
    plt.show()
"""
heat_map("expression_matrix_addRowand_Column_all0to5years_output.csv")
#heat_map("NDDCandidatesSubNetwork-12_output_ITC_averaged.csv")
#heat_map("ITC_STC_V1C_CB_PFC_0to5year_NDDs_correlation_map.csv")