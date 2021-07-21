import itertools

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

plt.rcParams.update({'font.size': 6})


def convert_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset = dataset.dropna(axis='index')  # drop rows which contain missing values
    dataset = dataset.drop(labels=['inchikey', 'cluster', 'split', 'fp', 'label'], axis='columns').copy()
    dataset = dataset.set_index('smiles')

    # print('Warning: dropping rows with positive scores')
    # positive_score_row = (dataset > 0).any(axis=1)
    # return dataset.loc[~positive_score_row]

    print('Warning: setting positive scores to 0.0')
    dataset[dataset > 0.0] = 0.0
    return dataset


def main() -> None:
    path = '/home/gregor/projects/cambridge/research/dockgym/downloads/big_dataset.tsv'
    dataset = convert_dataset(pd.read_csv(path, sep='\t'))

    pearson_rs = {(a, b): stats.pearsonr(dataset[a], dataset[b])[0]
                  for a, b in itertools.combinations(dataset.columns, 2)}

    with pd.option_context('display.max_rows', None):
        print(pd.Series(pearson_rs, name='pearson_r').sort_values())

    pairs = [('ADRB1', 'ADRB2'), ('MAOB', 'CYP3A4')]
    offset = 0.75

    fig_size = 2.5
    num_cols = len(pairs)
    fig, axes = plt.subplots(nrows=1, ncols=num_cols, figsize=(num_cols * fig_size, fig_size), constrained_layout=True)

    for ax, (target_a, target_b) in zip(axes, pairs):
        lower = np.min(dataset[[target_a, target_b]].to_numpy())
        upper = np.max(dataset[[target_a, target_b]].to_numpy())

        hb = ax.hexbin(
            dataset[target_a],
            dataset[target_b],
            gridsize=(15, 15),
            extent=(lower - 2 * offset, upper + 2 * offset, lower - 2 * offset, upper + 2 * offset),
            bins='log',
            cmap='Blues',
            linewidths=0.1,
        )

        ax.plot(
            (lower - 2 * offset, upper + 2 * offset),
            (lower - 2 * offset, upper + 2 * offset),
            linestyle='dashed',
            zorder=1,
            color='black',
            alpha=0.5,
        )

        ax.set_xlim(lower - offset, upper + offset)
        ax.set_ylim(lower - offset, upper + offset)

        ax.set_xlabel(target_a)
        ax.set_ylabel(target_b)

        fig.colorbar(hb, ax=ax)
        # cb.set_label('counts')

        print(f'Correlation coefficient {target_a}, {target_b}: {pearson_rs[(target_a, target_b)]}')

    fig.savefig('correlation.pdf')


if __name__ == '__main__':
    main()
