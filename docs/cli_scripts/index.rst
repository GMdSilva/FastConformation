Command-Line Interface (CLI) Guide
==================================

FastEnsemble provides command-line tools for performing various tasks such as MSA generation, ensemble prediction, and analysis. This section provides more detailed documentation for each CLI script.

.. toctree::
   :maxdepth: 2
   :caption: CLI Tools

   MSA_Usage
   Predictions_Usage
   Analysis_Usage

Overview
--------

Below, you'll find documentation for each tool, including the command-line arguments, usage examples, and best practices.

Available CLI Tools
-------------------

- **MSA Generation** (`msa.rst`):
  
  The `msa` script handles Multiple Sequence Alignment (MSA) generation using tools like `jackhmmer` and `mmseqs2`.

- **Ensemble Prediction** (`predict_ensemble.rst`):
  
  The `predict_ensemble` script is used for making ensemble predictions based on the generated MSAs.

- **Analysis** (`analysis.rst`):
  
  The `analysis` script provides tools for analyzing the predicted structures, including RMSD and TM-score calculations, PCA clustering, and more.

