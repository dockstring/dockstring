import argparse
from typing import Sequence

import pandas as pd
import sklearn.metrics as metrics


def parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--balanced', help='path to balanced TSV file', required=True)
    parser.add_argument('--unbalanced', help='path to unbalanced TSV file', required=True)
    return parser.parse_args(args=args)


def enrichment_factor(labels: Sequence[bool], scores: Sequence[float], top_k: int) -> float:
    data = pd.DataFrame({'label': labels, 'score': scores}).sort_values(by='score', ascending=False)
    num_actives_before = int(sum(data['label']))
    num_actives_after = int(sum(data['label'][:top_k]))
    rate_before = num_actives_before / data.shape[0]
    rate_after = num_actives_after / top_k
    return rate_after / rate_before


def prepare_dataset(dataset: pd.DataFrame) -> pd.DataFrame:
    dataset['label'] = dataset['label'] == 'A'  # note: some entries don't have a label -> false
    dataset = dataset.dropna(axis=0).copy()  # discard rows without score
    dataset['score'] = dataset['score'] * -1  # for the docking score, higher is better
    return dataset


def main():
    args = parse_args()

    # Balanced: AUC
    balanced_set = prepare_dataset(pd.read_csv(args.balanced, sep='\t'))
    # Unbalanced: AP, EF
    unbalanced_set = prepare_dataset(pd.read_csv(args.unbalanced, sep='\t'))

    metrics_list = []
    for prop in ['score', 'logp', 'qed']:
        auc = balanced_set.groupby('target').apply(lambda g: metrics.roc_auc_score(g['label'], g[prop]))
        ap = unbalanced_set.groupby('target').apply(lambda g: metrics.average_precision_score(g['label'], g[prop]))
        ef = unbalanced_set.groupby('target').apply(
            lambda g: enrichment_factor(labels=g['label'], scores=g[prop], top_k=200))

        score_metrics = pd.DataFrame({'auc': auc, 'ap': ap, 'ef': ef}).sort_values(by='auc', ascending=False)
        score_metrics['prop'] = prop
        metrics_list.append(score_metrics)

    quality_metrics = pd.concat(metrics_list)
    quality_metrics = quality_metrics.reset_index().set_index(['target', 'prop']).sort_index()
    quality_metrics.to_csv('quality_metrics_aug.tsv', sep='\t', header=True, index=True)

    quality_metrics.to_latex('table.tex', float_format='%.2f', sparsify=True, longtable=True)


if __name__ == '__main__':
    main()
