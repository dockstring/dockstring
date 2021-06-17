from typing import Sequence

import pandas as pd
import sklearn.metrics as metrics

from dockstring.utils import get_resources_dir


def enrichment_factor(labels: Sequence[bool], scores: Sequence[float], top_k: int) -> float:
    data = pd.DataFrame({'label': labels, 'score': scores}).sort_values(by='score', ascending=False)
    num_actives_before = int(sum(data['label']))
    num_actives_after = int(sum(data['label'][:top_k]))
    rate_before = num_actives_before / data.shape[0]
    rate_after = num_actives_after / top_k
    return rate_after / rate_before


def prepare_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset['label'] = dataset['label'] == 'A'  # note: some entries don't have a label -> false
    dataset = dataset[~dataset['score'].isna()].copy()  # discard entries without score
    dataset['score'] = dataset['score'] * -1  # higher is better
    return dataset


def main():
    # Balanced: AUC
    balanced_path = get_resources_dir() / 'data' / 'balanced_scores_1000_actives_1000_inactives.tsv'
    balanced_set = prepare_dataset(pd.read_csv(balanced_path, sep='\t'))

    auc = balanced_set.groupby('target').apply(lambda g: metrics.roc_auc_score(g['label'], g['score']))

    # Unbalanced: AP, EF
    unbalanced_path = get_resources_dir() / 'data' / 'unbalanced_scores_200_actives_3800_inactives.tsv'
    unbalanced_set = prepare_dataset(pd.read_csv(unbalanced_path, sep='\t'))

    ap = unbalanced_set.groupby('target').apply(lambda g: metrics.average_precision_score(g['label'], g['score']))
    ef = unbalanced_set.groupby('target').apply(
        lambda g: enrichment_factor(labels=g['label'], scores=g['score'], top_k=200))

    quality_metrics = pd.DataFrame({'auc': auc, 'ap': ap, 'ef': ef}).sort_values(by='auc', ascending=False)

    quality_metrics.to_csv(get_resources_dir() / 'data' / 'quality_metrics.tsv', sep='\t', header=True, index=True)


if __name__ == '__main__':
    main()
