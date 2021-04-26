import tempfile
from pathlib import Path
import sys
import os
import platform

def get_path_to_lib():
    return Path(__file__).parent.parent / 'lib'

def load_target(name):
    return Target(name)

class Target():
    def __init__(self,name):
        self._name = name
        # Create temporary directory where the PDB, PDBQT and conf files for the target will be saved
        self._tmp_dir_handle = tempfile.TemporaryDirectory()
        self._tmp_dir = Path(self._tmp_dir_handle.name).resolve()
        # Copy PDB, PDBQT and conf files to the temporary directory
        # TODO Check if PDB is needed or not ANSWER Maybe yes for visualization in pymol? Since PDBQT is not visualized the same way
        # TODO Create PDBQT for receptors
        pdb_reference = Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'receptors' / (self.name + '_receptor.pdb')
        shutil.copyfile(pdb_reference,self._pdb)
        pdbqt_reference = Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'receptors' / (self.name + '_receptor.pdbqt')
        shutil.copyfile(pdbqt_reference,self._pdbqt)
        conf_reference = Path(sys.modules[self.__class__.__module__].__file__).parent.parent / 'receptors' / (self.name + '_conf.txt')
        shutil.copyfile(conf_reference,self._conf)


    def __del__(self):
        self._tmp_dir_handle.cleanup()

    @property
    def name(self):
        return self._name
    # TODO
    @property
    def _pdb(self):
        return Path(self._tmp_dir / (self.name + '_receptor.pdb'))

    # TODO
    @property
    def _pdbqt(self):
        return Path(self._tmp_dir / (self.name + '_receptor.pdbqt'))

    @property
    def _conf(self):
        return Path(self._tmp_dir / (self.name + '_conf.txt'))


    def dock(self,mol=None,seed=95390476):
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





