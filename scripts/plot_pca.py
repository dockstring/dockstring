from typing import Optional

import pandas as pd
from matplotlib import pyplot as plt
from rdkit.Chem import Descriptors, AllChem
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


def convert_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset = dataset.dropna(axis=0)  # discard columns without score
    dataset = dataset.drop(labels=['inchikey', 'cluster', 'split', 'fp', 'label'], axis='columns').copy()
    dataset[dataset.select_dtypes(include=['number']).columns] *= -1
    return dataset


def parse_smiles(smiles: str) -> Optional[AllChem.Mol]:
    mol = AllChem.MolFromSmiles(smiles)
    if not mol:
        print(f'Cannot parse {smiles}')
    return mol


def compute_logp(mol: Optional[AllChem.Mol]) -> Optional[float]:
    return Descriptors.MolLogP(mol)


def compute_qed(mol: Optional[AllChem.Mol]) -> Optional[float]:
    return AllChem.QED.qed(mol)


def main() -> None:
    path = '/home/gregor/projects/cambridge/research/dockgym/downloads/big_dataset.tsv'
    dataset = convert_dataset(pd.read_csv(path, sep='\t', low_memory=False))
    dataset = dataset.loc[:10_000]
    print(dataset.columns)

    # Add logp and qed
    dataset['mol'] = dataset['smiles'].apply(parse_smiles)
    dataset['logp'] = dataset['mol'].apply(compute_logp)
    dataset['qed'] = dataset['mol'].apply(compute_qed)
    dataset = dataset.drop(labels=['smiles', 'mol'], axis='columns')
    print(dataset)

    # Compute PCA
    array = dataset.to_numpy(dtype=float).T  # [n_smiles, n_targets]
    print(array.shape)

    pca = PCA(n_components=2)
    pca.fit(array)

    print(pca)
    projected = pca.transform(array)

    # Plot
    fig_width = 5.0
    fig_height = 3.0
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(fig_width, fig_height), constrained_layout=True)

    ax.scatter(projected[:, 0], projected[:, 1])

    fig.show()


if __name__ == '__main__':
    main()
