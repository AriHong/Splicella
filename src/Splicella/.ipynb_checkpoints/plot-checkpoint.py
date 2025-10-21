import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def plot_geneisonum(spcDCI, ax):
    if ax is None:
        fig, ax = plt.subplots(1,1)
    sns.histplot(spcDCI.gene_iso_num[1], binrange=(.5, 5.5), bins=5,
            color='black', ax=ax)
    ax.set_xticks(range(1, 6))
    ax.set_xlabel('#isoform per gene')
    ax.set_ylabel('#gene')

    return None


