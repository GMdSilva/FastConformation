====================================
MSA Generation Tool Documentation
====================================

This documentation provides an overview of two Python scripts designed to automate the process of generating multiple sequence alignments (MSA) using either `jackhmmer` or `mmseqs2`.

.. contents::
   :local:
   :depth: 2

=======================================
JackHMMeR Multiple Sequence Alignment
=======================================

Overview
--------

This script is designed to generate a Multiple Sequence Alignment (MSA) using `jackhmmer`. It processes a target sequence provided in FASTA format and outputs the resulting MSA in a format compatible with tools like `colabfold_batch`.

**Key Features:**
- Configurable via command-line arguments or a JSON configuration file.
- Supports the use of a RAM disk for improved performance.
- Saves the resulting MSA in a specified output directory.

Usage
-----

To run the script, use the following command:

.. code-block:: bash

   jackhmmer_msa --config_file <path_to_config> --sequence_path <path_to_fasta> --output_path <output_dir> [optional arguments]

**Command-Line Arguments:**
- `--config_file`: Path to the JSON configuration file.
- `--sequence_path`: Path to the FASTA file containing the target sequence.
- `--output_path`: Path to save the results.
- `--jobname`: The job name (optional).
- `--homooligomers`: Number of copies of the protein (optional).
- `--use_ramdisk`: Whether to use a RAM disk for the process (optional, requires root access).

Example
-------

.. code-block:: bash

   jackhmmer_msa --config_file config.json --sequence_path input.fasta --output_path ./results

=====================================
mmseqs2 Multiple Sequence Alignment
=====================================

Overview
--------

This script is designed to generate a Multiple Sequence Alignment (MSA) using `mmseqs2`. Similar to the `jackhmmer` script, it processes a target sequence provided in FASTA format and outputs the resulting MSA.

**Key Features:**
- Configurable via command-line arguments or a JSON configuration file.
- Saves the resulting MSA and associated files in the specified output directory.

Usage
-----

To run the script, use the following command:

.. code-block:: bash

   mmseqs2_msa --config_file <path_to_config> --sequence_path <path_to_fasta> --output_path <output_dir> [optional arguments]

**Command-Line Arguments:**
- `--config_file`: Path to the JSON configuration file.
- `--sequence_path`: Path to the FASTA file containing the target sequence.
- `--output_path`: Path to save the results.
- `--jobname`: The job name (optional).

Example
-------

.. code-block:: bash

   mmseqs2_msa --config_file config.json --sequence_path input.fasta --output_path ./results

=====================
Configuration Files
=====================

Both scripts allow configuration through a JSON file. Below is an example configuration file:

.. code-block:: json

   {
       "sequence_path": "input.fasta",
       "output_path": "./results",
       "jobname": "prediction_run",
       "homooligomers": 1,
       "use_ramdisk": false,
       "tmp_dir": "./tmp"
   }

====================
Output
====================

Both scripts save the resulting MSA in the specified output directory. The MSA file can be used for downstream analyses, such as structure prediction or further alignment refinement.
