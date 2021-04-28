from guatask import Task, run_task
import sys
from dockgym import load_target
import dockgym.src.utils as utils
from rdkit import Chem

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
        mol = target.dock('CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)OP(=O)(O)OC(C)CNC(=O)CCC4(C(C5C6(C(C(C(=N6)C(=C7C(C(C(=N7)C=C8C(C(C(=N8)C(=C4[N-]5)C)CCC(=O)N)(C)C)CCC(=O)N)(C)CC(=O)N)C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+2]')
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
        mol = target.dock('CN1CCc2cccc3c2[C@H]1Cc1ccc(OCCCNC(=O)CCCCCCCCCCCCCCCCCCC(=O)NCCCOc2ccc4c(c2O)-c2cccc5c2[C@@H](C4)N(C)CC5)c(O)c1-3')
        if mol is None:
            print('Fail')
        else:
            print('Pass')
        __import__('pdb').set_trace()
        print('hi')


    def load_output(self):
        raise NotImplementedError('Task.load_output not implemented')


if __name__ == '__main__':
    task = sys.argv[1]
    #__import__('pdb').set_trace()
    if task in globals():
        run_task(globals()[task])
    else:
        raise RuntimeError('Task {} is not defined.'.format(task))
