import itertools

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

plt.rcParams.update({'font.size': 6})

colors = [
    '#1f77b4',  # muted blue
    '#d62728',  # brick red
    '#ff7f0e',  # safety orange
    '#2ca02c',  # cooked asparagus green
    '#9467bd',  # muted purple
    '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf',  # blue-teal
]


def convert_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset = dataset.dropna(axis='index')  # drop rows which contain missing values
    dataset = dataset.drop(labels=['inchikey', 'cluster', 'split', 'fp', 'label'], axis='columns').copy()
    dataset = dataset.set_index('smiles')

    print('Warning: dropping rows with positive scores')
    positive_score_row = (dataset > 0).any(axis=1)
    return dataset.loc[~positive_score_row]


def main() -> None:
    path = '/home/gregor/projects/cambridge/research/dockgym/downloads/big_dataset.tsv'
    dataset = convert_dataset(pd.read_csv(path, sep='\t'))

    pearson_rs = {(a, b): stats.pearsonr(dataset[a], dataset[b])[0]
                  for a, b in itertools.combinations(dataset.columns, 2)}

    with pd.option_context('display.max_rows', None):
        print(pd.Series(pearson_rs, name='pearson_r').sort_values())

    pairs = [('ADRB1', 'ADRB2'), ('MAOB', 'CYP3A4')]

    fig_size = 2.5
    num_cols = len(pairs)
    fig, axes = plt.subplots(nrows=1, ncols=num_cols, figsize=(num_cols * fig_size, fig_size), constrained_layout=True)

    for ax, (target_a, target_b) in zip(axes, pairs):
        ax.scatter(dataset[target_a], dataset[target_b], color='black', s=7, alpha=0.75)
        print(f'Correlation coefficient {target_a}, {target_b}: {pearson_rs[(target_a, target_b)]}')

        ax.set_xlabel(target_a)
        ax.set_ylabel(target_b)

        # Plot XY line
        limits = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
        ax.plot(limits, limits, linestyle='dashed', zorder=0, color=colors[0])

    fig.savefig('correlation.pdf')


if __name__ == '__main__':
    main()
