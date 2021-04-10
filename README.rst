====================
Nextflow-lcep
====================

This pipeline serves as a reference implementation for machine learning models build with `mlf-core <mlf-core.com>`_ using the mlflow-xgboost template.
The herefore trained mlf-core model can be found `here <https://github.com/mlf-core/lcep>`_ together with the corresponding `python package <https://github.com/mlf-core/lcep-package>`_.

Requirements
~~~~~~~~~~~~~

All predictions performed during pipeline execution are run on the GPU. To facilitate GPU usage you require CUDA and the NVIDIA container toolkit installed.

Usage
~~~~~~~

Build the Docker container using ``docker build -t nextflow/lcep:1.0.0 .``.

Next, run the pipeline via
``nextflow run main.nf -with-docker --pathways data/pathways_ensembl_cancer_2021.csv``
