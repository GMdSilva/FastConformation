Ensemble Prediction
=====================

``predict_ensemble`` allows to predict different protein conformations starting from an input MSA by run AlphaFold2 using the ColabFold implementation with different subsampling parameters.
Below is a detailed description of each argument and how to use them effectively.

Command-Line Arguments
-----------------------

- **--config_file** (str):
  
  Path to the configuration file. If not provided, the script will default to `config.json` in the current directory.
  
- **--jobname** (str):
  
  The name of the job. This is used to organize output directories and files.

- **--msa_path** (str):
  
  Path to the `.a3m` MSA file. If not provided, the script will automatically generate it based on the `output_path` and `jobname`.

- **--output_path** (str):
  
  Directory path where the prediction results will be saved.

- **--seq_pairs** (str):
  
  A list of `[max_seq, extra_seq]` pairs in the format `[[max_seq1, extra_seq1], [max_seq2, extra_seq2], ...]`. This defines the sequence pairing strategy for the predictions.

- **--seeds** (int, nargs='+'):
  
  Specifies the number of predictions to run. The default is 10.

- **--save_all** (bool):
  
  Outputs a pickled files of all the output.

- **--platform** (str):
  
  The platform to run the predictions on, either `cpu` or `gpu`. The default is `cpu`.

- **--subset_msa_to** (int):
  
  Subset the input MSA to the specified number of sequences.

- **--msa_from** (str):
  
  The MSA building tool used to generate the input MSA. Available options are `jackhmmer` or `mmseqs2`.

Usage Examples
--------------

Example 1: Using a Configuration File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have a configuration file named `config.json`, you can run the script as follows:

.. code-block:: bash

  predict_ensemble --config_file config.json
