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

    align
    make_fastq
    make_bwa_index


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
    utils.merge_fastq_files