from guatask import Task, run_task
import sys
from dockgym import load_target

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
        target.view()
    def load_output(self):
        raise NotImplementedError('Task.load_output not implemented')


if __name__ == '__main__':
    task = sys.argv[1]
    #__import__('pdb').set_trace()
    if task in globals():
        run_task(globals()[task])
    else:
        raise RuntimeError('Task {} is not defined.'.format(task))
