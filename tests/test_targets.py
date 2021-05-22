import math
import os
import tempfile
from pathlib import Path

import pytest
import rdkit.Chem as Chem

from dockstring import list_all_target_names, load_target, DockingError
from dockstring.utils import (smiles_to_mol, write_embedded_mol_to_pdb, embed_mol, read_mol_from_pdb, protonate_pdb,
                              check_vina_output, parse_scores_from_output, canonicalize_smiles)


class TestLoader:
    def test_load_all_targets(self):
        names = list_all_target_names()
        assert len(names) == 58
        assert all(isinstance(name, str) for name in names)
        assert all(load_target(name) for name in names)

    def test_wrong_target(self):
        with pytest.raises(DockingError):
            load_target('does_not_exist')


lysine_smiles = 'C(CCN)C[C@@H](C(=O)O)N'
aspartic_acid_smiles = 'C([C@@H](C(=O)O)N)C(=O)O'
alanine_smiles = 'CC(C(=O)O)N'


class TestConversions:
    def test_convert_string_success(self):
        assert smiles_to_mol('C')

    def test_convert_string_fail(self):
        with pytest.raises(DockingError):
            smiles_to_mol('not_a_mol')

    def test_charged_mol(self):
        with pytest.raises(DockingError):
            smiles_to_mol('CCC(=O)O{-1}')

    def test_write_fail(self):
        mol = smiles_to_mol(lysine_smiles)
        with tempfile.NamedTemporaryFile(suffix='.pdb') as f:
            with pytest.raises(DockingError):
                write_embedded_mol_to_pdb(mol, ligand_pdb=f.name)

    @pytest.mark.parametrize('smiles', [
        lysine_smiles,
        'O=C1N(C=2N=C(OC)N=CC2N=C1C=3C=CC=CC3)C4CC4',
    ])
    def test_read_write_pdb(self, smiles):
        mol = smiles_to_mol(smiles)
        embedded_mol = embed_mol(mol, seed=1)

        # Added Hs
        assert Chem.MolToSmiles(embedded_mol) != Chem.MolToSmiles(mol)

        with tempfile.NamedTemporaryFile(suffix='.pdb') as f:
            write_embedded_mol_to_pdb(embedded_mol, ligand_pdb=f.name)
            read_mol = read_mol_from_pdb(f.name)

        # Hs are gone
        assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(read_mol)

    @pytest.mark.parametrize('smiles,charge_ph7', [
        (lysine_smiles, 1),
        (alanine_smiles, 0),
        (aspartic_acid_smiles, -1),
    ])
    def test_protonation(self, smiles, charge_ph7):
        mol = smiles_to_mol(smiles)
        embedded_mol = embed_mol(mol, seed=1)

        with tempfile.NamedTemporaryFile(suffix='.pdb') as pdb_file:
            write_embedded_mol_to_pdb(embedded_mol, ligand_pdb=pdb_file.name)
            protonate_pdb(pdb_file.name)
            read_mol = read_mol_from_pdb(pdb_file.name)

        charges = tuple(sum(atom.GetFormalCharge() for atom in m.GetAtoms()) for m in (mol, read_mol))
        assert charges == (0, charge_ph7)


resources_dir = Path(os.path.dirname(os.path.realpath(__file__))) / 'resources'


class TestParser:
    def test_score_parser(self):
        vina_output = resources_dir / 'vina.out'
        assert check_vina_output(vina_output) is None

        scores = parse_scores_from_output(vina_output)
        expected = [-4.7, -4.6, -4.5, -4.5, -4.4, -4.4, -4.4, -4.3, -4.3]
        assert len(scores) == len(expected)
        assert all(math.isclose(a, b) for a, b in zip(scores, expected))


class TestEmbedding:
    @pytest.mark.parametrize('smiles', [
        'S(C=1C=2NCC(CC2C=C(C1)C(=O)OCC)(C)C)(N[C@H](C(N3CCC(CC3)CCF)=O)CC4=NC=5C=CC=CC5S4)(=O)=O',
        'N1(C([C@@H](C2=CSC=C2)NC(OC(C)(C)C)=O)=O)[C@H](C(N[C@H](C=O)CCCN=C(N)N)=O)CCC1',
        'P(OC1=CC=C(C[C@@H](C(N[C@H](C(=O)NCCC=2C=CC=CC2)CCC(O)=O)=O)NC(=O)C)C=C1)(=O)(O)O',
        r'OC1=C(C=CC(=CC/C=C(/CC/C=C(/CCC=C(C)C)\C)\C)C)C(=C(O)C=2C1=CC=CC2)C',
    ])
    def test_difficult(self, smiles: str):
        canonical_smiles = canonicalize_smiles(smiles)
        mol = smiles_to_mol(canonical_smiles)
        embedded_mol = embed_mol(mol, seed=1)
        assert embedded_mol.GetNumConformers() == 1

    @pytest.mark.parametrize('smiles', [
        'OC=1N(C(O)=C2C1C3=C4C(=C2CC3O)C=CC=C4)CC5=CC=CC=C5',
    ])
    def test_impossible(self, smiles: str):
        canonical_smiles = canonicalize_smiles(smiles)
        mol = smiles_to_mol(canonical_smiles)
        with pytest.raises(DockingError):
            embed_mol(mol, seed=1)


class TestRefinement:
    def test_successful_refinement(self):
        pass

    def test_no_ff_parameters(self):
        pass


class TestDocking:
    def test_simple_docking(self):
        target = load_target('ABL1')

        smiles_1 = 'CCO'
        energy_1, _ = target.dock(smiles_1)
        assert math.isclose(energy_1, -2.4)

        smiles_2 = 'CC'
        energy_2, _ = target.dock(smiles_2)
        assert math.isclose(energy_2, -1.8)

    # Test different SMILES representations of lysine
    @pytest.mark.parametrize('smiles', [lysine_smiles, 'NCCCC[C@H](N)C(=O)O'])
    def test_charged(self, smiles):
        target = load_target('CYP3A4')
        energy, aux = target.dock(smiles)
        assert math.isclose(energy, -4.7)

        charge = sum(atom.GetFormalCharge() for atom in aux['ligand'].GetAtoms())
        assert charge == 1

    def test_pdbqt_to_pdb_error(self):
        target = load_target('CYP3A4')
        score, aux = target.dock('O=C1N(C=2N=C(OC)N=CC2N=C1C=3C=CC=CC3)C4CC4')
        scores = [-9.1, -8.5, -8.2, -8.2, -8.0, -7.9, -7.8, -7.8, -7.7]
        assert aux['ligand'].GetNumConformers() == 9
        assert len(aux['scores']) == len(scores)
        assert all(math.isclose(a, b) for a, b in zip(aux['scores'], scores))

    def test_multiple_molecules(self):
        target = load_target('ABL1')
        with pytest.raises(DockingError):
            target.dock('C.C')
        with pytest.raises(DockingError):
            target.dock('C.CO')
