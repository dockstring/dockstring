from dockstring.utils import get_resources_dir
import pandas as pd
import sklearn.metrics as metrics

balanced_set = pd.read_csv(get_resources_dir() / 'data' / 'balanced_scores_1000_actives_1000_inactives.tsv', sep='\t')
unbalanced_set = pd.read_csv(get_resources_dir() / 'data' / 'unbalanced_scores_200_actives_3800_inactives.tsv', sep='\t')

def enrichment_factor(labels,scores,top_k):
    data = pd.DataFrame([labels,scores]).transpose().sort_values(by='score',ascending=False)
    num_actives_before = int(sum(data['label']))
    num_actives_after = int(sum(data['label'].iloc[:top_k]))
    rate_before = num_actives_before / data.shape[0]
    rate_after = num_actives_after / top_k
    ef = rate_after / rate_before
    return ef

targets = balanced_set['target'].unique()
quality_metrics = pd.DataFrame(index=targets,columns=['auc','ap','ef'])

for each_target in targets:
    # Load balanced and unbalanced sets
    # - make integer labels
    # - make scores positive (so that higher is better)
    # - remove missing values
    each_balanced_set = balanced_set.loc[balanced_set['target'] == each_target]
    each_balanced_set['label'] = (each_balanced_set['label'] == 'A').astype(int)
    each_balanced_set['score'] = each_balanced_set['score'] * -1
    each_balanced_set = each_balanced_set.loc[~each_balanced_set['score'].isnull()]
    each_unbalanced_set = unbalanced_set.loc[unbalanced_set['target'] == each_target]
    each_unbalanced_set['label'] = (each_unbalanced_set['label'] == 'A').astype(int)
    each_unbalanced_set['score'] = each_unbalanced_set['score'] * -1
    each_unbalanced_set = each_unbalanced_set.loc[~each_unbalanced_set['score'].isnull()]
    # Compute area under a ROC curve and average precision
    labels = each_balanced_set['label']
    scores = each_balanced_set['score']
    auc = metrics.roc_auc_score(labels,scores)
    # Compute enrichment factor
    labels = each_unbalanced_set['label']
    scores = each_unbalanced_set['score']
    ap = metrics.average_precision_score(labels,scores)
    ef = enrichment_factor(labels=labels,scores=scores,top_k=200)
    # Save to dataframe
    quality_metrics.loc[each_target,['auc','ap','ef']] = [auc,ap,ef]


quality_metrics = quality_metrics.sort_values(by='auc',ascending=False)
quality_metrics.to_csv(get_resources_dir() / 'data' / 'quality_metrics.tsv',sep='\t',header=True,index=True)
