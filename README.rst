**circconvert**
======================================================================

**circRNA database conversion tool**

.. figure:: https://raw.githubusercontent.com/dieterich-lab/circtools/master/docs/img/circtools_200px.png
   :alt: circconvert

|downloads| |pypi|

Introduction
-------------

Circular RNAs (circRNAs) originate through back-splicing events from linear primary transcripts, are resistant to exonucleases, typically not polyadenylated, and have been shown to be highly specific for cell type and developmental stage. Although few circular RNA molecules have been shown to exhibit miRNA sponge function, for the vast majority of circRNAs however, their function is yet to be determined.

The prediction of circular RNAs is a multi-stage bioinformatics process starting with raw sequencing data and usually ending with a list of potential circRNA candidates which, depending on tissue and condition may contain hundreds to thousands of potential circRNAs. While there already exist a number of tools for the prediction process (e.g. `DCC <https://github.com/dieterich-lab/DCC>`__ and `CircTest <https://github.com/dieterich-lab/CircTest>`__), publicly available downstream analysis tools are rare.

We developed **circtools**, a modular, Python3-based framework for circRNA-related tools that unifies several functionalities in single command line driven software. The command line follows the `circtools subcommand` standard that is employed in samtools or bedtools. Currently, circtools includes modules for detecting and reconstructing circRNAs,
a quick check of circRNA mapping results, RBP enrichment screenings, circRNA primer design, statistical testing, and an exon usage module.

Installation
------------

The ``circtools`` package is written in Python3 (>=3.4), two modules, namely ``detect`` and ``reconstruct`` also require a working Python 2 installation (>=2.7). It requires only a small number of external dependencies, namely standard bioinformatics tools:

Installation is managed through ``python3 setup.py install``. No sudo
access is required if the installation is executed with ``--user`` which
will install the package in a user-writeable folder. The binaries should
be installed to ``/home/$user/.local/bin/`` in case of Debian-based
systems.

``circtools`` was developed and tested on Debian Jessie but should also
run with any distribution.

The installation requires running python on the command line:

::

    git clone https://github.com/jakobilab/circconvert.git
    cd circconvert
    python3 setup.py install --verbose --user


Modules
-------

Circconvert currently offers two modules:

convert
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



query
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. |downloads| image:: https://pepy.tech/badge/circtools
    :alt: Python Package Index Downloads
    :scale: 100%
    :target: https://pepy.tech/project/circtools

.. |pypi| image:: https://badge.fury.io/py/circtools.svg
    :alt: Python package version
    :scale: 100%
    :target: https://badge.fury.io/py/circtools
