import math
import os
import tempfile
from pathlib import Path

import pytest
from rdkit.Chem import AllChem as Chem

from dockstring import list_all_target_names, load_target
from dockstring.errors import (
    DockstringError,
    EmbeddingError,
    FormatConversionError,
    ParsingError,
    PoseProcessingError,
    SanityError,
    StructureOptimizationError,
)
from dockstring.utils import (
    canonicalize_smiles,
    check_vina_output,
    embed_mol,
    parse_affinities_from_output,
    protonate_mol,
    refine_mol_with_ff,
    smiles_to_mol,
    write_mol_to_mol_file,
)


class TestLoader:
    def test_load_all_targets(self):
        names = list_all_target_names()
        assert len(names) == 58
        assert all(isinstance(name, str) for name in names)
        assert all(load_target(name) for name in names)

    def test_wrong_target(self):
        with pytest.raises(DockstringError):
            load_target('does_not_exist')


lysine_smiles = 'C(CCN)C[C@@H](C(=O)O)N'
aspartic_acid_smiles = 'C([C@@H](C(=O)O)N)C(=O)O'
alanine_smiles = 'CC(C(=O)O)N'


class TestConversions:
    def test_convert_string_success(self):
        assert smiles_to_mol('C')

    def test_convert_string_fail(self):
        with pytest.raises(ParsingError):
            smiles_to_mol('not_a_mol')

    def test_charged_mol(self):
        smiles_to_mol('CCC(=O)[O-]')
        smiles_to_mol('CC(C)(C)CC(C)(C)C1=CC=C(C=C1)OCCOCC[N+](C)(C)CC2=CC=CC=C2')

    def test_write_fail(self):
        mol = smiles_to_mol(lysine_smiles)
        with tempfile.NamedTemporaryFile(suffix='.mol') as f:
            with pytest.raises(DockstringError):
                write_mol_to_mol_file(mol, mol_file=f.name)

    @pytest.mark.parametrize('smiles,charge_ph7', [
        (lysine_smiles, 1),
        (alanine_smiles, 0),
        (aspartic_acid_smiles, -1),
    ])
    def test_protonation(self, smiles, charge_ph7):
        mol = smiles_to_mol(smiles)
        protonated_mol = protonate_mol(mol, pH=7.4)
        charges = tuple(sum(atom.GetFormalCharge() for atom in m.GetAtoms()) for m in (mol, protonated_mol))
        assert charges == (0, charge_ph7)


resources_dir = Path(os.path.dirname(os.path.realpath(__file__))) / 'resources'


class TestParser:
    def test_affinities_parser(self):
        vina_output = resources_dir / 'vina.out'
        check_vina_output(vina_output)

        affinities = parse_affinities_from_output(vina_output)
        expected = [-4.7, -4.6, -4.5, -4.5, -4.4, -4.4, -4.4, -4.3, -4.3]
        assert len(affinities) == len(expected)
        assert all(math.isclose(a, b) for a, b in zip(affinities, expected))


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
        with pytest.raises(EmbeddingError):
            embed_mol(mol, seed=1)


class TestRefinement:
    @pytest.mark.parametrize('smiles', [
        'C(NC=1C=CC=CC1)(C2=CC=C(C=C2)C3=CC=C(C(NC=4C=CC=CC4)=O)C=C3)=O',
        'C=1(C=CC(=CC1)NC(=O)C2CCN(CC2)C(=O)C3=CC=C(C=C3)F)C=4C=CC(=CC4)NC(=O)C5CCN(CC5)C(C6=CC=C(C=C6)F)=O',
    ])
    def test_successful_refinement(self, smiles):
        canonical_smiles = canonicalize_smiles(smiles)
        mol = smiles_to_mol(canonical_smiles)
        embedded_mol = embed_mol(mol, seed=1)

        props = Chem.MMFFGetMoleculeProperties(embedded_mol)
        initial_energy = Chem.MMFFGetMoleculeForceField(embedded_mol, props).CalcEnergy()
        refined_mol = refine_mol_with_ff(embedded_mol)
        final_energy = Chem.MMFFGetMoleculeForceField(refined_mol, props).CalcEnergy()

        assert final_energy < initial_energy

    @pytest.mark.parametrize('smiles', [
        'S(O)(O)(N=C(C(N1CCN(CC1)C)C=2SC=CC2)C)=CC',
        'S(N1[C@@H](C(O)=NO)CC(C1)=NS(O)(=O)C)(C2=CC=C(C=C2)OC)(=O)=O',
    ])
    def test_no_ff_parameters(self, smiles):
        canonical_smiles = canonicalize_smiles(smiles)
        mol = smiles_to_mol(canonical_smiles)
        embedded_mol = embed_mol(mol, seed=1)
        with pytest.raises(StructureOptimizationError):
            refine_mol_with_ff(embedded_mol)

    @pytest.mark.parametrize('smiles', [
        'S=1(O)(O)=CC=2C(=NN(C2NC(=O)C=3OC=4C(C3)=CC=CC4)C5=CC=C(C=C5)C)C1',
        'ClC1=CC=C(N2N=C3C(C=S(O)(O)=C3)=C2NC(=O)C4CC4)C=C1',
        'S=1(O)(O)=CC2=C(N(N=C2C1)C3=CC=C(F)C=C3)NC(=O)C45CC6CC(C4)CC(C5)C6',
        'S=1(O)(O)=CC=2C(=NN(C2NC(=O)C=3OC=4C(C3)=CC=CC4)C5=CC=C(C=C5)C)C1',
    ])
    def test_kekulization(self, smiles):
        # With UFF these molecules can be optimized, not with MMFF though (raises KekulizationError)
        canonical_smiles = canonicalize_smiles(smiles)
        mol = smiles_to_mol(canonical_smiles)
        embedded_mol = embed_mol(mol, seed=1)
        assert refine_mol_with_ff(embedded_mol)


def test_simple_docking():
    """Simple docking test, set to run in all cases."""
    target = load_target('ABL1')

    smiles_1 = 'CCO'
    energy_1, _ = target.dock(smiles_1)
    assert energy_1 is not None and math.isclose(energy_1, -2.4)

    smiles_2 = 'CC'
    energy_2, _ = target.dock(smiles_2)
    assert energy_2 is not None and math.isclose(energy_2, -1.8)


@pytest.mark.slow
class TestDocking_Slow:
    """Run a large number of slower docking tests."""

    # Test different SMILES representations of lysine
    @pytest.mark.parametrize('smiles', [lysine_smiles, 'NCCCC[C@H](N)C(=O)O'])
    def test_charged(self, smiles):
        target = load_target('CYP3A4')
        energy, aux = target.dock(smiles)
        assert energy is not None and math.isclose(energy, -4.6)

        charge = sum(atom.GetFormalCharge() for atom in aux['ligand'].GetAtoms())
        assert charge == 1

    @pytest.mark.parametrize('seed, score', [(0, -4.6), (1, -4.6), (2, -4.7)])
    def test_seeds(self, seed: int, score: float):
        target = load_target('CYP3A4')
        energy, aux = target.dock(lysine_smiles, seed=seed)
        assert energy is not None and math.isclose(energy, score)

        charge = sum(atom.GetFormalCharge() for atom in aux['ligand'].GetAtoms())
        assert charge == 1

    @pytest.mark.parametrize(
        'smiles, charge, energy',
        [
            ('[H][N+]1=CC=CC=C1', 0, -4.2),  # pyridinium
            ('N1=CC=CC=C1', 0, -4.2),  # pyridine
            ('C[N+]1=CC=CC=C1', 1, -4.3),
            ('CC(O)=O', -1, -3.0),  # acetic acid
            ('CC([O-])=O', -1, -3.0),  # acetic acid
        ])
    def test_different_charges(self, smiles, charge, energy):
        target = load_target('CYP3A4')
        docking_energy, aux = target.dock(smiles)
        total_charge = sum(atom.GetFormalCharge() for atom in aux['ligand'].GetAtoms())
        assert total_charge == charge
        assert docking_energy is not None and math.isclose(energy, docking_energy)

    def test_positive_score(self):
        target = load_target('AR')
        smiles = r'O1/C(=N\CCC=2C=3C(NC2)=CC=CC3)/[C@]4(N(C(=O)C5C=6[C@H]4[C@H]7[C@@H](CC6[C@@H]8[C@H]([C@H]5C)C(=O)N(C8=O)C9=CC=CC=C9)C(=O)N(C7=O)C%10=CC=CC=C%10)C1=O)CC%11=CC=CC=C%11'
        score, aux = target.dock(smiles)
        assert score is not None and score >= 0.0

    def test_mol_to_pdbqt_error(self):
        # works for any target actually
        target = load_target('CYP3A4')
        with pytest.raises(FormatConversionError):
            target.dock(
                'S(=O)(=O)(N1C=C(C2(C[C@H](NC(OC(C)(C)C)=O)C(OC)=O)C3=C(NC2=O)C=CC=C3)C=4C1=CC=CC4)CC[Si](C)(C)C')

    def test_pdbqt_to_pdb_error(self):
        target = load_target('CYP3A4')
        score, aux = target.dock('O=C1N(C=2N=C(OC)N=CC2N=C1C=3C=CC=CC3)C4CC4')
        affinities = [-9.0, -8.8, -8.4, -8.3, -8.3, -7.9, -7.7, -7.7, -7.6]
        assert aux['ligand'].GetNumConformers() == 9
        assert len(aux['affinities']) == len(affinities)
        assert all(math.isclose(a, b) for a, b in zip(aux['affinities'], affinities))

    def test_chiral_centers(self):
        target = load_target('CYP3A4')

        score, aux = target.dock('[H]C(C)(F)Cl')
        assert score is not None and math.isclose(score, -3.1)
        assert Chem.MolToSmiles(aux['ligand']) == 'C[C@H](F)Cl'

        score, aux = target.dock('C[C@H](F)Cl')
        assert score is not None and math.isclose(score, -3.1)
        assert Chem.MolToSmiles(aux['ligand']) == 'C[C@H](F)Cl'

    def test_bond_stereo(self):
        target = load_target('CYP3A4')

        smiles = r'C\C=C\C'  # E
        score, aux = target.dock(smiles)
        assert score is not None and math.isclose(score, -3.4)
        assert aux['ligand'].GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREOE

        smiles = r'C/C=C\C'  # Z
        score, aux = target.dock(smiles)
        assert score is not None and math.isclose(score, -3.5)
        assert aux['ligand'].GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREOZ

        smiles = 'CC=CC'  # unspecified
        score, aux = target.dock(smiles)
        assert score is not None and math.isclose(score, -3.5)
        assert aux['ligand'].GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREOZ

    @pytest.mark.parametrize('target_name, ligand', [
        ('MAPK1', 'C1=CC=C2C3=C(NC2=C1)[C@H](N4C(=O)CN(C(=O)[C@H]4C3)C)C5=CC=C6OCOC6=C5'),
        ('MAOB', 'C1(=CC(=C(C=C1)N2CCOCC2)F)N3C[C@H](CNC(C)=O)OC3=O'),
        ('ABL1', 'S(=O)(=O)(N(CC)CC)C1=C(SC)C=CC(=C1)C'),
    ])
    def test_additional_chiral_ligands(self, target_name: str, ligand: str):
        target = load_target(target_name)
        assert target.dock(ligand)

    def test_multiple_molecules(self):
        target = load_target('ABL1')
        with pytest.raises(SanityError):
            target.dock('C.C')
        with pytest.raises(SanityError):
            target.dock('C.CO')

    def test_radicals(self):
        target = load_target('ABL1')
        with pytest.raises(SanityError):
            target.dock('C[CH2]')

        with pytest.raises(SanityError):
            target.dock('C[CH]')

    @pytest.mark.parametrize(
        'target_name, ligand_smiles',
        [
            # ('ABL1', 'BrC12CC3(CC(C1)CC(C3)C2)CC(=O)NCC4=CC=CC=C4'),  # too many bonds
            ('ABL1', 'BrC1=CC(P(OCC)(OCC)=O)(NS(=O)(=O)C2=CC=CC=C2)C3=C(C1=O)C=CC=C3'),  # multiple fragments
        ])
    def test_bond_assignment_fails(self, target_name: str, ligand_smiles: str):
        target = load_target('ABL1')
        with pytest.raises(PoseProcessingError):
            target.dock(ligand_smiles)

    @pytest.mark.parametrize('target_name, ligand', [
        ('ABL1',
         'S(=O)(=O)(N(C[C@@H]1OCCCC[C@@H](OC=2C(C(=O)N(C[C@H]1C)[C@@H](CO)C)=CC(NC(=O)NC=3C=CC(F)=CC3)=CC2)C)C)C=4SC=CC4'
         ),
    ])
    def test_atom_valence_exception(self, target_name: str, ligand: str):
        target = load_target(target_name)
        assert target.dock(ligand)
