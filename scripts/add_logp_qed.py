import argparse
from typing import Optional

import pandas as pd
from rdkit.Chem import Descriptors, AllChem


def parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', help='path to input TSV file', required=True)
    return parser.parse_args(args=args)


def parse_smiles(smiles: str) -> Optional[AllChem.Mol]:
    mol = AllChem.MolFromSmiles(smiles)
    if not mol:
        print(f'Cannot parse {smiles}')
    return mol


def compute_logp(mol: Optional[AllChem.Mol]) -> Optional[float]:
    if mol:
        return Descriptors.MolLogP(mol)
    return None


def compute_qed(mol: Optional[AllChem.Mol]) -> Optional[float]:
    if mol:
        return AllChem.QED.qed(mol)
    return None


def main():
    args = parse_args()

    print(f'Reading file: {args.input}')
    df = pd.read_csv(args.input, sep='\t')

    df['mols'] = df['smiles'].apply(parse_smiles)
    df['logp'] = df['mols'].apply(compute_logp)
    df['qed'] = df['mols'].apply(compute_qed)

    df = df.drop('mols', axis='columns')

    output_file_name = 'augmented_dataset.tsv'
    print(f'Writing file: {output_file_name}')
    df.to_csv(output_file_name, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
