Ensemble Prediction
====================

Command-Line Arguments
-----------------------

predict_ensemble can be executed with the following command-line arguments:

- `--config_file` (str): Path to the configuration file. If not provided, the script will look for `config.json` in the current directory.
  
- `--jobname` (str): The name of the job. This will be used to organize output directories and files.
  
- `--msa_path` (str): Path to the `.a3m` MSA file. If not provided, it will be auto-generated based on the `output_path` and `jobname`.
  
- `--output_path` (str): Directory path where the prediction results will be saved.
  
- `--seq_pairs` (str): List of `[max_seq, extra_seq]` pairs in the format `[[max_seq1, extra_seq1], [max_seq2, extra_seq2], ...]`.
  
- `--seeds` (int, nargs='+'): Number of predictions to run. The default is 10.
  
- `--save_all` (bool): Flag to save all results. If not provided, only the most relevant results will be saved.
  
- `--platform` (str): Platform to run the predictions (`cpu` or `gpu`). The default is `cpu`.
  
- `--subset_msa_to` (int): Subsets the input MSA to the specified number of sequences. Useful when running out of RAM due to deep MSAs.
  
- `--msa_from` (str): MSA building tool used to build the input MSA. Choices are `jackhmmer` or `mmseqs2`.

Usage Examples
--------------

Example 1: Using a configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have a configuration file named `config.json`, you can simply run the script as follows:

```bash
python -m fast_ensemble.predict_ensemble --config_file config.json