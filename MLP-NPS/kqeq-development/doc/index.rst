kQEq
================================

kQEq is a Python implementation of the kernel Charge Equilibration method (more `here <https://doi.org/xxx>`_). This allows training and using physics-based machine-learning models for predicting charge distributions in molecules and materials.

This documentation contains a brief summary of the theory, a tutorial and the API documentation.

Installation
============

This package has been tested with python 3.9. External dependencies are `numpy <https://numpy.org>`_ (for linear algebra), `ase <https://wiki.fysik.dtu.dk/ase/>`_ (for handling structural data) and `DScribe <https://singroup.github.io/dscribe/latest/>`_ (for calculating atomic environment representations). The kQEq package and its dependencies can be installed using ``pip install .`` 


Contents
========
.. toctree::

    theory
    tutorial
    API


  
Changelog
=========
 - 0.1:
    - Initial release.
