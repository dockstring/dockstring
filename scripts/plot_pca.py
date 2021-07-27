import argparse

import pandas as pd
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA

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


def parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', help='path to input TSV file', required=True)
    return parser.parse_args(args=args)


def convert_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset = dataset.dropna(axis=0)  # discard columns without score
    dataset = dataset.drop(labels=['inchikey', 'smiles', 'cluster', 'split', 'fp', 'label'], axis='columns').copy()
    return dataset


def main() -> None:
    args = parse_args()
    print(f'Reading file: {args.dataset}')
    dataset = convert_dataset(pd.read_csv(args.dataset, sep='\t', low_memory=False))

    array = dataset.to_numpy(dtype=float).T  # [n_smiles, n_targets]

    # Standardize data
    array -= array.mean(axis=0)
    array /= array.std(axis=0)

    # Compute PCA
    pca = PCA(n_components=2)
    pca.fit(array)

    projected = pca.transform(array)

    # Plot
    fig_size = 5.50107 / 2  # inches, NeurIPS template
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_size, fig_size), constrained_layout=True)

    ax.scatter(projected[:, 0], projected[:, 1], s=7, color='black')

    # for i, label in enumerate(dataset.columns):
    #     plt.annotate(label, (projected[i, 0], projected[i, 1]))

    ax.tick_params(
        axis='both',
        which='both',
        left=False,
        top=False,
        right=False,
        bottom=False,
        labelbottom=False,
        labelleft=False,
    )
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')

    fig.savefig('pca.pdf')


if __name__ == '__main__':
    main()
