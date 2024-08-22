Usage
=====

This section provides detailed instructions on how to use FastEnsemble, including how to run the GUI and various CLI tools provided by the package.

Installing on Linux
===================

To install FastEnsemble on a Linux system, follow these steps:

1. Open a terminal.
2. Run the installation script by typing:

   .. code-block:: bash

      bash ./install.sh

   Running the installation script will open the GUI.

3. To run the GUI again after exiting, type:

   .. code-block:: bash

      poetry run run_gui

   Or you can also run it as:

   .. code-block:: bash

      fast_ensemble/fast_ensemble-conda/bin/python -m fast_ensemble.gui.run_gui

Installing on Windows (WSL)
===========================

To run FastEnsemble on Windows using WSL, ensure that WSL and Miniconda are installed. You can install WSL by following these steps:

1. Open PowerShell, Windows Command Prompt, or Terminal in administrator mode (Right-click > Run as administrator).
2. Type the following command to install WSL:

   .. code-block:: bash

      wsl --install

3. Restart your machine.

Running FastEnsemble on Windows
===============================

1. Download the `install.sh` file from the GitHub repository.
2. Open the terminal and type `WSL` to enter the WSL environment.
3. Navigate to the WSL file system by typing:

   .. code-block:: bash

      cd ~/

4. Move the downloaded `install.sh` file to this directory:

   .. code-block:: bash

      mv /location/of/file/install.sh ~/

   Replace `/location/of/file/` with the path to where your `install.sh` file is located.

5. Run the installation script by typing:

   .. code-block:: bash

      bash ./install.sh

   Running the installation script will open the GUI.

6. To run the GUI again after exiting, type:

   .. code-block:: bash

      fast_ensemble/fast_ensemble-conda/bin/python -m fast_ensemble.gui.run_gui

   Note that you must be in the WSL file system (`~/`) to run this command successfully.

Running the GUI
===============

Once FastEnsemble is installed, you can access its graphical user interface (GUI) by running the following command:

.. code-block:: bash

   $run_gui

Alternatively, you can start the GUI using:

.. code-block:: bash

   fast_ensemble/fast_ensemble-conda/bin/python -m fast_ensemble.gui.run_gui

This GUI allows you to perform MSA generation, AF2 prediction, and analysis to predict different protein conformations using MSA subsampling.

Using the Command-Line Interface (CLI)
======================================

FastEnsemble provides several CLI tools for different tasks, including MSA generation, prediction, and analysis. Below is a list of available commands with brief descriptions.

**MSA Generation:**

- **jackhmmer_msa**: Generate MSA using `jackhmmer`.

  .. code-block:: bash

     $jackhmmer_msa --config_file <path_to_config> --sequence_path <path_to_fasta> --output_path <output_dir>


- **mmseqs2_msa**: Generate MSA using `mmseqs2`.

  .. code-block:: bash

     $mmseqs2_msa --config_file <path_to_config> --sequence_path <path_to_fasta> --output_path <output_dir>


**Prediction:**

- **predict_ensemble**: Run ensemble predictions.

  .. code-block:: bash

     $predict_ensemble --config_file <path_to_config>

- **decaf_init**: Initialize DECAF predictions.

  .. code-block:: bash

     $decaf_init --config_file <path_to_config>

**Analysis:**

- **rmsd_mode1d**: Analyze RMSD in 1D mode.

  .. code-block:: bash

     $rmsd_mode1d --config_file <path_to_config>

- **rmsd_mode2d**: Analyze RMSD in 2D mode.

  .. code-block:: bash

     $rmsd_mode2d --config_file <path_to_config>

- **tmscore_mode1d**: Analyze TM-score in 1D mode.

  .. code-block:: bash

     $tmscore_mode1d --config_file <path_to_config>

- **tmscore_mode2d**: Analyze TM-score in 2D mode.

  .. code-block:: bash

     $tmscore_mode2d --config_file <path_to_config>

- **pca_clustering**: Perform PCA clustering on the predicted structures.

  .. code-block:: bash

     $pca_clustering --config_file <path_to_config>

- **rmsf_plddt**: Calculate RMSF and pLDDT for the predicted structures.

  .. code-block:: bash

     $rmsf_plddt --config_file <path_to_config>

- **save_traj**: Save trajectories from the analysis.

  .. code-block:: bash

     $save_traj --config_file <path_to_config>

For more detailed instructions on how to use each tool, refer to the respective CLI documentation sections provided in this guide.

Additional Notes
================

- Ensure that all necessary configuration files are correctly set up as described in their respective sections.
- For more detailed instructions and examples, visit the specific CLI documentation sections.
