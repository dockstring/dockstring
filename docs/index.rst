.. dockstring documentation master file, created by
   sphinx-quickstart on Tue Sep 28 16:04:37 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dockstring's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Dockstring is very easy to use.

.. code-block:: python
   :linenos:

   from dockstring import load_target

   # Choose which protein you want to dock
   target = load_target("DRD2")

   # Specify molecule as a smiles string
   smiles = "CCCC=O"

   # Run the docking simulation
   score, info = target.dock(smiles)
   print(score)

That's it! For more details see our API below:

.. autoclass:: dockstring.target.Target
    :members:
    :show-inheritance:

.. automodule:: dockstring.errors
    :members:
    :show-inheritance:

.. automodule:: dockstring
    :members:
    :undoc-members:
    :show-inheritance:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
