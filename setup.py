import os
from typing import Dict

from setuptools import setup, find_packages


def readme() -> str:
    with open('README.md') as f:
        return f.read()


version_dict = {}  # type: Dict[str, str]
with open(os.path.join('dockstring', '_version.py')) as file:
    exec(file.read(), version_dict)

setup(
    name='dockstring',
    version=version_dict['__version__'],
    description='Package for easy ligand preparation and docking with AutoDock Vina against '
    'a highly-curated selection of 58 protein targets. One of the three '
    'components of the DOCKSTRING bundle (package, dataset, and benchmark tasks).',
    long_description=readme(),
    classifiers=['Programming Language :: Python :: 3.7'],
    python_requires='>=3.7',
    packages=find_packages(),
    include_package_data=True,
    package_data={'dockstring': ['resources/targets/*', 'resources/bin/*']},
    install_requires=[],
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest'],
)
