import tempfile
from pathlib import Path
import sys
import os
import shutil
import platform
from rdkit.Chem import AllChem as Chem
from .utils import robust, DockingError, get_num_conf

def load_target(name):
    return Target(name)

class Target():
    def __init__(self,name,random_seed=974528263):
        self._name = name
        self.random_seed = random_seed
        # Create temporary directory where the PDB, PDBQT and conf files for the target will be saved
        self._tmp_dir_handle = tempfile.TemporaryDirectory()
        # Copy PDB, PDBQT and conf files to the temporary directory
        # TODO Create PDBQT for receptors
        pdb_reference = self._receptors_dir  / (self.name + '_receptor.pdb')
        shutil.copyfile(pdb_reference,self._pdb)
        pdbqt_reference = self._receptors_dir / (self.name + '_receptor.pdbqt')
        shutil.copyfile(pdbqt_reference,self._pdbqt)
        conf_reference = self._receptors_dir / (self.name + '_conf.txt')
        shutil.copyfile(conf_reference,self._conf)


    def __del__(self):
        self._tmp_dir_handle.cleanup()

    @property
    def name(self):
        return self._name
    # Paths to important locations as property methods
    @property
    def _tmp_dir(self):
        return Path(self._tmp_dir_handle.name).resolve()
    @property
    def _receptors_dir(self):
        return Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'receptors'
    @property
    def _bin_dir(self):
        return Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'bin'
    @property
    def _bin_dir(self):
        return Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'bin'
    @property
    def _pdb(self):
        return Path(self._tmp_dir / (self.name + '_receptor.pdb'))
    @property
    def _pdbqt(self):
        return Path(self._tmp_dir / (self.name + '_receptor.pdbqt'))
    @property
    def _conf(self):
        return Path(self._tmp_dir / (self.name + '_conf.txt'))

    # Submethods for docking a ligand
    def _smiles_2_mol(self,smiles):
        mol =  Chem.MolFromSmiles(smiles)
        if mol is None:
            raise DockingError(f'Docking of molecule  {smiles}  failed during the ' \
                    'conversion of SMILES to RDKit molecule object.')
        return mol

    def _inchi_2_mol(self,inchi):
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            raise DockingError(f'Docking of molecule  {inchi}  failed during the ' \
                    'conversion of SMILES to RDKit molecule object.')
        return mol

    def _mol_2_embedding(self,mol):
        # Will attempt to find 3D coordinates 10 times with different random seeds
        n_attempts = 10
        # Set first random seed for reproducibility
        random_seed = self._dock_random_seed
        for i in range(n_attempts):
            # Simple approach to get subsequent random seeds from first random seed
            random_seed = (random_seed - (-1)**(i % 2) * i * random_seed) % 766523564
            # Always add hydrogens in order to get a sensible 3D structure, and remove them later
            mol = Chem.AddHs(mol)
            Chem.EmbedMolecule(mol,randomSeed=random_seed)
            mol = Chem.RemoveHs(mol)
            # If at least one conformation has been obtained, don't try more and break out of the loop.
            # Otherwise, keep trying to hopefully generate a ligand conformation
            if get_num_conf(mol) > 0:
                break
        # If not a single conformation is obtained in all the attempts, raise an error and return None.
        # Otherwise, return the molecule with the conformation
        if get_num_conf(mol) == 0:
            raise DockingError(f'Docking of molecule  {self._mol_id}  failed during the ' \
                                'ligand conformation generation.')
        else:
            return mol
        '''
        - mol: either a SMILES string, an inchikey or a RDKit molecule object
        - seed: integer random seed for reproducibility
        '''
        # Docking with Vina is performed in a temporary directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            # Define necessary filenames
            pdb_filename   = (tmp_dir / 'ligand.pdb').resolve()
            pdbqt_filename = (tmp_dir / 'ligand.pdbqt').resolve()
            vina_logfile   = (tmp_dir / 'ligand.log').resolve()
            vina_outfile   = (tmp_dir / 'ligand.out').resolve()
            # Define paths to dependencies
            lib = get_path_to_lib()
            prep_ligand = lib / 'prepare_ligand4.py'
            system = platform.system()
            if system == 'Linux':
                vina_binary = lib / 'vina_linux'
                pythonsh_binary = lib / 'pythonsh' # TODO Check if pythonsh is platform-dependent, or it is the same for Linux, Mac, Windows, etc
            elif system == 'Windows':
                pass
            elif system == 'Darwin':
                pass
            # TODO Save all the temporary files in this directory
            # Using tempfile makes the temporary file creation portable for all systems
            print(tmp_dir)

        # TODO Put the vina and pythonsh executables for Mac and Windows under "lib". And also the prepare_ligand4.py
        # TODO Check whether the current system is Linux, Mac or Windows with platform.system()
        #      https://stackoverflow.com/questions/22321397/python-os-name-return-nt-on-windows-7
        # TODO Put all the receptors and configuration files under "receptors"
        # TODO Put all the calculated scores (and maybe the poses too?) under "data". What should be the format?
        # - Plain text for the scores and smiles?
        # - What format for the poses?

    def info(self):
        '''
        Print some info about the target. Maybe just a string?
        '''
        pass

    def view(self, search_box=True):
        '''
        Start pymol and view the receptor and the search box.
        '''
        with open(self._conf,'r') as f:
            # Extract the search box information
            for line in f:
                if   'center_x' in line:
                    center_x = float(line.split()[2])
                elif 'center_y' in line:
                    center_y = float(line.split()[2])
                elif 'center_z' in line:
                    center_z = float(line.split()[2])
                elif 'size_x'   in line:
                    size_x   = float(line.split()[2])
                elif 'size_y'   in line:
                    size_y   = float(line.split()[2])
                elif 'size_z'   in line:
                    size_z   = float(line.split()[2])

            receptor = str(self._pdb.resolve())
            # TODO Remove ligand from the main view() method, make it a part of the docked ligands' method
            ligand = f'/home/mgarort/repos/dockgym/playground/crystal_ligands/{self.name}/crystal_ligand.mol2'
            pymol_view_search_box = Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'utils' / 'view_search_box.py'
            command = f'pymol -R {pymol_view_search_box} {receptor} {ligand}'
            if search_box:
                command += f" -d 'view_search_box center_x={center_x}, center_y={center_y}, center_z={center_z}, \
                              size_x={size_x}, size_y={size_y}, size_z={size_z}'"
            command += ' > /dev/null'
            # TODO Change to subprocess
            os.system(command)





