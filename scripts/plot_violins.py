import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from dockstring.utils import get_resources_dir


def prepare_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset = dataset.dropna(axis=0)  # discard rows without score
    return dataset[dataset['score'] < 0.0]


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

# yapf: disable
targets_1 = [
    'ABL1', 'ADAM17', 'ADRB1', 'ADRB2', 'AKT2', 'MAOB', 'CASP3', 'DHFR', 'ESR2', 'PTK2', 'FGFR1', 'HMGCR',
    'HSP90AA1', 'KIT', 'MAPKAPK2', 'MAP2K1', 'NOS1', 'PARP1', 'PDE5A', 'PPARD', 'PGR', 'PTPN1', 'ROCK1',
    'AKT1', 'AR', 'CDK2', 'CSF1R', 'ESR1',
]
targets_2 = [
    'NR3C1', 'IGF1R', 'JAK2', 'LCK', 'MET', 'MMP13', 'PTGS2', 'PPARA', 'PPARG', 'REN', 'ADORA2A', 'ACHE', 'BACE1',
    'CA2', 'CYP2C9', 'CYP3A4', 'HSD11B1', 'DPP4', 'DRD2', 'DRD3', 'EGFR', 'F10', 'GBA', 'MAPK1', 'MAPK14', 'PLK1',
    'SRC', 'THRB', 'F2', 'KDR'
]
# yapf: enable


def main() -> None:
    path = get_resources_dir() / 'data' / 'unbalanced_scores_200_actives_3800_inactives.tsv'
    dataset = prepare_dataset(pd.read_csv(path, sep='\t'))

    assert len(set(targets_1 + targets_2)) == len(set(dataset['target'].unique()))

    fig_width = 5.0
    fig_height = 3.0

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(fig_width, fig_height), constrained_layout=True)

    for ax, targets in zip(axes, (targets_1, targets_2)):
        data = [sorted(dataset.loc[dataset['target'] == target, 'score']) for target in targets]
        parts = ax.violinplot(data, showextrema=False, showmedians=True)

        for body in parts['bodies']:
            body.set_facecolor(colors[0])
            body.set_edgecolor('black')

        parts['cmedians'].set_edgecolor(colors[0])

        ax.set_xticks(np.arange(1, len(targets) + 1))
        ax.set_xticklabels(targets, rotation=90)

    fig.savefig('violin.pdf')


if __name__ == '__main__':
    main()
