import math
import tempfile

import pytest
import rdkit.Chem as Chem

from dockstring import list_all_target_names, load_target, DockingError
from dockstring.utils import (smiles_or_inchi_to_mol, write_embedded_mol_to_pdb, embed_mol, read_mol_from_pdb,
                              protonate_pdb, convert_pdb_to_pdbqt, convert_pdbqt_to_pdb)


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


class TestConversions:
    def test_convert_string_success(self):
        assert smiles_or_inchi_to_mol('C')

    def test_convert_string_fail(self):
        with pytest.raises(DockingError):
            smiles_or_inchi_to_mol('not_a_mol')

    def test_write_fail(self):
        mol = smiles_or_inchi_to_mol(lysine_smiles)
        with tempfile.NamedTemporaryFile(suffix='.pdb') as f:
            with pytest.raises(DockingError):
                write_embedded_mol_to_pdb(mol, ligand_pdb=f.name)

    @pytest.mark.parametrize('smiles', [
        lysine_smiles,
        'O=C1N(C=2N=C(OC)N=CC2N=C1C=3C=CC=CC3)C4CC4',
    ])
    def test_read_write_pdb(self, smiles):
        mol = smiles_or_inchi_to_mol(smiles)
        embedded_mol = embed_mol(mol, seed=1)

        # Added Hs
        assert Chem.MolToSmiles(embedded_mol) != Chem.MolToSmiles(mol)

        with tempfile.NamedTemporaryFile(suffix='.pdb') as f:
            write_embedded_mol_to_pdb(embedded_mol, ligand_pdb=f.name)
            read_mol = read_mol_from_pdb(f.name)

        # Hs are gone
        assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(read_mol)

    def test_protonation(self):
        mol = smiles_or_inchi_to_mol(lysine_smiles)
        embedded_mol = embed_mol(mol, seed=1)

        with tempfile.NamedTemporaryFile(suffix='.pdb') as pdb_file:
            write_embedded_mol_to_pdb(embedded_mol, ligand_pdb=pdb_file.name)
            protonate_pdb(pdb_file.name)

            # Check that conversion to PDBQT doesn't break anything
            with tempfile.NamedTemporaryFile(suffix='.pdbqt') as pdbqt_file:
                convert_pdb_to_pdbqt(pdb_file.name, pdbqt_file.name)
                convert_pdbqt_to_pdb(pdbqt_file.name, pdb_file.name)

            read_mol = read_mol_from_pdb(pdb_file.name)

        charges = tuple(sum(atom.GetFormalCharge() for atom in m.GetAtoms()) for m in (mol, read_mol))
        assert charges == (0, 2)


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

    def test_charged(self):
        target = load_target('CYP3A4')
        energy, aux = target.dock(lysine_smiles)
        assert math.isclose(energy, -4.6)

        charge = sum(atom.GetFormalCharge() for atom in aux['ligands'].GetAtoms())
        assert charge == 2

    def test_pdbqt_to_pdb_error(self):
        target = load_target('CYP3A4')
        result = target.dock('O=C1N(C=2N=C(OC)N=CC2N=C1C=3C=CC=CC3)C4CC4')
        # assert math.isclose(score, -9.1)  # cannot read PDBQT file properly
        assert result == (None, None)
