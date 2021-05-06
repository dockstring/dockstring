from guatask import Task, run_task
import sys
from dockgym import load_target
import dockgym.src.utils as utils
from rdkit import Chem
import pandas as pd


class Debug(Task):
    requires = []
    params = {}
    directory = 'debug_directory'
    subdirectory = 'debug_subdirectory'
    output_filename = 'NONE'
    debug = True

    def run(self):
        receptor = sys.argv[2]
        target = load_target(receptor)
        # mol = target.dock('C1=CC=CC=C1')
        # mol = target.dock('CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)OP(=O)(O)OC(C)CNC(=O)CCC4(C(C5C6(C(C(C(=N6)C(=C7C(C(C(=N7)C=C8C(C(C(=N8)C(=C4[N-]5)C)CCC(=O)N)(C)C)CCC(=O)N)(C)CC(=O)N)C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+2]')
        target.view()

    def load_output(self):
        raise NotImplementedError('Task.load_output not implemented')


class FailedEmbedProcessing(Task):
    requires = []
    params = {}
    directory = 'debug_directory'
    subdirectory = 'debug_subdirectory'
    output_filename = 'NONE'
    debug = True

    def run(self):
        receptor = sys.argv[2]
        target = load_target(receptor)
        print('Should print: Pass, Pass, Fail')
        # Pass
        mol = target.dock('C1=CC=CC=C1')
        if mol is None:
            print('Fail')
        else:
            print('Pass')
        mol = target.dock('CC1CC2C(C2(C)C)C3C=C(C(C4(C1(C3=O)C=C(C4O)C)O)O)CO')
        if mol is None:
            print('Fail')
        else:
            print('Pass')
        mol = target.dock(
            'CN1CCc2cccc3c2[C@H]1Cc1ccc(OCCCNC(=O)CCCCCCCCCCCCCCCCCCC(=O)NCCCOc2ccc4c(c2O)-c2cccc5c2[C@@H](C4)N(C)CC5)c(O)c1-3'
        )
        if mol is None:
            print('Fail')
        else:
            print('Pass')
        __import__('pdb').set_trace()
        print('hi')

    def load_output(self):
        raise NotImplementedError('Task.load_output not implemented')


class TestDocking(Task):
    requires = []
    params = {}
    directory = 'debug_directory'
    subdirectory = 'debug_subdirectory'
    output_filename = 'NONE'
    debug = True

    def run(self):
        receptor = sys.argv[2]
        target = load_target(receptor)
        # mol = target.dock('C1=CC=CC=C1', seed=28,logfile='hello',verbose=True)
        # mol = target.dock('C1=CC=CC=C1OBEF*@&#BR', seed=28,logfile='hello',verbose=True)
        mol = target.dock(
            'InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)',
            seed=28,
            logfile='hello',
            verbose=True)
        # score, other = target.dock('CC=C(C)C(=O)OC1C(=CC23C1(C(C(=CC(C2=O)C4C(C4(C)C)CC3C)CO)O)O)C',seed=10, num_cpu=1)
        # target.view()
    def load_output(self):
        raise NotImplementedError('Task.load_output not implemented')


class TestAllPDBQTs(Task):
    requires = []
    # yapf: disable
    params = {'receptors': ['ABL1'     , 'ADAM17' , 'ADRB1'    , 'ADRB2'  , 'AKT2'   , 'MAOB'  ,
                            'CASP3'    , 'DHFR'   , 'ESR2'     , 'PTK2'   , 'FGFR1'  , 'HMGCR' ,
                            'HSP90AA1' , 'KIT'    , 'MAPKAPK2' , 'MAP2K1' , 'NOS1'   , 'PARP1' ,
                            'PDE5A'    , 'PPARD'  , 'PGR'      , 'PTPN1'  , 'ROCK1'  , 'AKT1'  ,
                            'AR'       , 'CDK2'   , 'CSF1R'    , 'ESR1'   , 'NR3C1'  , 'IGF1R' ,
                            'JAK2'     , 'LCK'    , 'MET'      , 'MMP13'  , 'PTGS2'  , 'PPARA' ,
                            'PPARG'    , 'REN'    , 'ADORA2A'  , 'ACHE'   , 'BACE1'  , 'CA2'   ,
                            'CYP2C9'   , 'CYP3A4' , 'HSD11B1'  , 'DPP4'   , 'DRD2'   , 'DRD3'  ,
                            'EGFR'     , 'F10'    , 'GBA'      , 'MAPK1'  , 'MAPK14' , 'PLK1'  ,
                            'SRC'      , 'THRB'   , 'F2'       , 'KDR']}
    # yapf: enable
    directory = 'tests'
    subdirectory = 'test_all_pdbqts'
    output_filename = 'scores.tsv'
    input_filename = '100_sample_smiles.txt'
    debug = True

    def run(self):
        # Initialize dataframe
        receptors = self.params['receptors']
        with open(self.input_filepath, 'r') as f:
            smiles = [line.strip().strip('\n') for line in f]
        scores = pd.DataFrame(index=smiles, columns=receptors)
        # Compute scores and save them in dataframe
        for each_receptor in receptors:
            target = load_target(each_receptor)
            for each_smiles in smiles:
                (score, other) = target.dock(each_smiles)
                scores.loc[each_smiles, each_receptor] = score
                print(each_receptor, each_smiles, score)
        # Save dataframe
        scores.to_csv(self.output_filepath, sep='\t')

    def load_output(self):
        pd.read_csv(self.output_filepath, sep='\t')


if __name__ == '__main__':
    task = sys.argv[1]
    #__import__('pdb').set_trace()
    if task in globals():
        run_task(globals()[task])
    else:
        raise RuntimeError('Task {} is not defined.'.format(task))
