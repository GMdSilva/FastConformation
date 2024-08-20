
Command-Line Interface (CLI) Usage
===================================

This module is used to perform TM-Score, RMSD (Root Mean Square Deviation) and RMSF (Root Mean Square Fluctuation) analyses on protein structure predictions. The analysis is configured using a JSON configuration file or command line arguments.

Running 2D TM-Score Analysis
----------------------------

The `run_2d_tmscore_analysis` function performs a 2D TM-Score analysis on the provided protein structure predictions. This analysis compares the structural similarity between pairs of sequences and reference structures.

**Example:**

.. code-block:: bash

   poetry run tmscore_mode2d --config_file config.json --output_path /path/to/output --predictions_path /path/to/predictions --jobname my_job

**Command Line Arguments:**

- `--config_file`: Path to the JSON configuration file (default: config.json).
- `--output_path`: Directory to save the results.
- `--predictions_path`: Directory containing the predictions.
- `--jobname`: Name of the job.
- `--seq_pairs`: A list of [max_seq, extra_seq] pairs used for predictions.
- `--starting_residue`: Starting residue index for reindexing (optional).
- `--slice_predictions`: Slice range of predictions to analyze (optional).
- `--engine`: The engine used to generate predictions (e.g., AlphaFold2, OpenFold).
- `--ref2d1`: First reference structure for TM-Score calculations (optional).
- `--ref2d2`: Second reference structure for TM-Score calculations (optional).
- `--n_stdevs`: Number of standard deviations to consider for point closeness (optional).
- `--n_clusters`: Number of clusters for TM-Score analysis (optional).

**Python Usage:**

You can also call the function directly in Python:

.. code-block:: python

   from fast_ensemble.ensemble_analysis.analysis_utils import load_config
   from script_name import run_2d_tmscore_analysis

   config = load_config('config.json')
   run_2d_tmscore_analysis(config)

Running 1D TM-Score Analysis
----------------------------

The `run_tmscore_analysis` function performs a 1D TM-Score analysis. This analysis identifies structural modes in the protein predictions and clusters them.

**Example:**

.. code-block:: bash

   poetry run tmscore_mode1d --config_file config.json --output_path /path/to/output --jobname my_job --engine AlphaFold2

**Command Line Arguments:**

- `--config_file`: Path to the JSON configuration file (default: config.json).
- `--output_path`: Directory to save the results.
- `--predictions_path`: Directory containing the predictions.
- `--jobname`: Name of the job.
- `--seq_pairs`: A list of [max_seq, extra_seq] pairs used for predictions.
- `--starting_residue`: Starting residue index for reindexing (optional).
- `--slice_predictions`: Slice range of predictions to analyze (optional).
- `--ref1`: Reference structure for TM-Score calculations (optional).
- `--engine`: The engine used to generate predictions (e.g., AlphaFold2, OpenFold).

**Python Usage:**

.. code-block:: python

   from fast_ensemble.ensemble_analysis.analysis_utils import load_config
   from script_name import run_tmscore_analysis

   config = load_config('config.json')
   run_tmscore_analysis(config)

Running 2D RMSD Analysis
------------------------

The `run_2d_rmsd_analysis` function performs a 2D RMSD analysis. This analysis measures the structural deviation between pairs of sequences and two reference structures in a 2D space.

**Example:**

.. code-block:: bash

   poetry run rmsd_mode2d --config_file config.json --output_path /path/to/output --predictions_path /path/to/predictions --jobname my_job

**Command Line Arguments:**

- `--config_file`: Path to the JSON configuration file (default: config.json).
- `--output_path`: Directory to save the results.
- `--mode_results`: Path to the mode results CSV file.
- `--jobname`: Name of the job.
- `--seq_pairs`: A list of [max_seq, extra_seq] pairs used for predictions.
- `--predictions_path`: Directory containing the predictions.
- `--engine`: The engine used to generate predictions (e.g., AlphaFold2, OpenFold).
- `--align_range`: Atom alignment range for RMSD calculations (optional).
- `--analysis_range`: Atom range for RMSD calculations after alignment (optional).
- `--analysis_range_name`: Name of the atom range (e.g., kinase core, helix 1, etc.).
- `--ref2d1`: First reference structure for RMSD calculations (optional).
- `--ref2d2`: Second reference structure for RMSD calculations (optional).
- `--n_stdevs`: Number of standard deviations to consider when calculating close points (optional).
- `--n_clusters`: Number of clusters to consider for RMSD analysis (optional).

**Python Usage:**

.. code-block:: python

   from fast_ensemble.ensemble_analysis.analysis_utils import load_config
   from script_name import run_2d_rmsd_analysis

   config = load_config('config.json')
   run_2d_rmsd_analysis(config)

Running 1D RMSD Analysis
------------------------

The `run_rmsd_analysis` function performs a 1D RMSD analysis. This analysis measures the structural deviation between sequences and a single reference structure.

**Example:**

.. code-block:: bash

   poetry run rmsd_mode1d --config_file config.json --output_path /path/to/output --jobname my_job --engine AlphaFold2

**Command Line Arguments:**

- `--config_file`: Path to the JSON configuration file (default: config.json).
- `--output_path`: Directory to save the results.
- `--predictions_path`: Directory containing the predictions.
- `--jobname`: Name of the job.
- `--seq_pairs`: A list of [max_seq, extra_seq] pairs used for predictions.
- `--starting_residue`: Starting residue index for reindexing (optional).
- `--align_range`: Atom alignment range for RMSF calculations (optional).
- `--analysis_range`: Atom range for RMSD calculations after alignment (optional).
- `--analysis_range_name`: Name of the atom range (e.g., kinase core, helix 1, etc.).
- `--ref1d`: Reference structure for RMSD calculations (optional).

**Python Usage:**

.. code-block:: python

   from fast_ensemble.ensemble_analysis.analysis_utils import load_config
   from script_name import run_rmsd_analysis

   config = load_config('config.json')
   run_rmsd_analysis(config)

Running RMSF Analysis
---------------------

The `run_rmsf_analysis` function performs RMSF analysis, which measures the flexibility of residues in the protein structure predictions.

**Example:**

.. code-block:: bash

   poetry run rmsf_plddt --config_file config.json --output_path /path/to/output --jobname my_job --engine AlphaFold2 --detect_mobile True

**Command Line Arguments:**

- `--config_file`: Path to the JSON configuration file (default: config.json).
- `--output_path`: Directory to save the results.
- `--predictions_path`: Directory containing the predictions.
- `--jobname`: Name of the job.
- `--seq_pairs`: A list of [max_seq, extra_seq] pairs used for predictions.
- `--engine`: The engine used to generate predictions (e.g., AlphaFold2, OpenFold).
- `--starting_residue`: Starting residue index for reindexing (optional).
- `--align_range`: Atom alignment range for RMSF calculations (optional).
- `--detect_mobile`: Boolean flag to detect mobile residue ranges (optional).
- `--peak_width`: RMSF peak width threshold for mobile residue range detection (optional).
- `--peak_prominence`: RMSF peak prominence threshold for mobile residue range detection (optional).
- `--peak_height`: RMSF peak height threshold for mobile residue range detection (optional).

**Python Usage:**

.. code-block:: python

   from fast_ensemble.ensemble_analysis.analysis_utils import load_config
   from script_name import run_rmsf_analysis

   config = load_config('config.json')
   run_rmsf_analysis(config)

Configuration
=============

The JSON configuration file should define the parameters necessary for each analysis. Here is an example configuration:

.. code-block:: json

   {
       "output_path": "/path/to/output",
       "predictions_path": "/path/to/predictions",
       "jobname": "my_job",
       "seq_pairs": [["seq1", "seq2"], ["seq3", "seq4"]],
       "engine": "AlphaFold2",
       "starting_residue": 1,
       "slice_predictions": "10:100",
       "align_range": "5-100",
       "detect_mobile": true,
       "peak_width": 3,
       "peak_prominence": 0.5,
       "peak_height": 1.0
   }

Each function will use the configuration parameters defined in the JSON file, but they can be overridden by command line arguments.
You can find example configs in the sample files.