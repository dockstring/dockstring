import tempfile
from pathlib import Path
import platform

def get_path_to_lib():
    return Path(__file__).parent.parent / 'lib'

class Target():
    def __init__(self):
        self._receptor = None # This should be a path to a temporary file with the pdbqt, which is copied from the repository

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

    def view(self):
        '''
        Start pymol and view the receptor and the search box.
        '''
        pass



