=============
API reference
=============

This page gives an overview of all public precellar objects, functions and
methods.

.. currentmodule:: precellar

Assay
~~~~~

.. autosummary::
    :toctree: _autosummary

    Assay

Core functions
~~~~~~~~~~~~~~

.. autosummary::
    :toctree: _autosummary

    make_genome_index
    align
    make_fragment
    make_fastq


Aligners
~~~~~~~~

.. autosummary::
    :toctree: _autosummary

    aligners.STAR
    aligners.BWAMEM2

Utilities
~~~~~~~~~

.. autosummary::
    :toctree: _autosummary

    utils.strip_barcode_from_fastq
    utils.bam_to_fastq