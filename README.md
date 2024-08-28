# FastEnsemble
![FastEnsemble Logo](background_logo.png)

FastEnsemble is a Python-based application that integrates MSA generation, structure prediction via AlphaFold 2 (AF2), and interactive analysis of protein conformational ensembles, all in one place. Uniquely, this tool enables researchers to leverage ML to generate protein conformations and analyze their populations without running MD simulations.

## Citation
FastEnsemble is based off of the research described in the manuscripts below.

Monteiro da Silva, G., Cui, J.Y., Dalgarno, D.C. et al. High-throughput prediction of protein conformational distributions with subsampled AlphaFold2. Nat Commun 15, 2464 (2024). https://doi.org/10.1038/s41467-024-46715-9 

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Download Sample Files](#download-sample-files)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Features

- **MSA Generation**: Automatically generates multiple sequence alignments (MSAs) from amino acid sequences using JACKHMMER and MMseqs2.
- **Structure Prediction**: Predicts protein structures using the ColabFold implementation of AlphaFold 2 (AF2).
- **Conformational Ensembles**: Generates alternative protein conformations through MSA subsampling.
- **Interactive Analysis**: Analyzes protein conformational ensembles and the effects of mutations on protein dynamics using our suite of analysis tools.
- **User-Friendly GUI**: Accessible through an intuitive graphical user interface suitable for non-programmers.

![alt text](image.png)

## Installation

To install FastEnsemble, download and run the installation script provided in the repository:

```bash
./install.sh
```

If you wish to install only the CLI version without the gui, run this script instead:

```bash
./install_cli.sh
```

This script will set up the necessary environment and dependencies.
At the end of the script we add the source path to the bashrc file so that you can use any of the commands without the need to activate the conda environment.

Next, close the terminal window and reopen a new one.

Now, run any of the commands reported in the Usage section.

## Dependencies

The installation script sets up the environment and installs all necessary dependencies. 

FastEnsemble relies on the following dependencies:

### Conda Packages

- `git`
- `python=3.10`
- `openmm==7.7.0`
- `pdbfixer`
- `kalign2=2.04`
- `hhsuite=3.3.0`
- `mmseqs2=15.6f452`
- `hmmer`
- `scikit-learn`
- `mdanalysis`
- `seaborn`
- `scipy`

### Python Packages (via pip)

- `PyQt5`
- `pandas`
- `pyqt`
- `matplotlib`
- `silence_tensorflow`
- `pyqtgraph`
- `colabfold`
- `silence_tensorflow`
- `pdb-tools`


## Documentation

Documentation can be found on [this ReadTheDocs page](https://fastensemble.readthedocs.io/en/main/index.html).


## OS Requirements

This package has been tested on Linux RedHat7 and Ubuntu 20.0.4.

## Usage
Before running any commands, download this file and run it:
```bash
./fast_ensemble.sh
```
### Running the GUI

To start the graphical user interface, execute the following command:

```bash
run_gui
```

This will launch the main application window, where you can access various functionalities such as submitting new jobs, checking job status, and viewing analysis logs.

### Running via Command-line

First, run 

```bash
fast_ensemble_init
```

Next, run any of the following, either specifying a config file path, or by specifying the parameters via command-line arguments. 

Sample config files and sample results are available via this link [Download Sample Files](https://drive.google.com/drive/folders/1ev5HfWVyMTBw3FRtKWxYaswuaIXvC1FS?usp=drive_link).

For more detailed instructions on how to use each tool, refer to the ReadTheDocs documentation.

**MSA Generation:**

- **jackhmmer_msa**: Generate MSA using `jackhmmer`.

    ```bash
    jackhmmer_msa --config_file <path_to_config>
    ```

- **mmseqs2_msa**: Generate MSA using `mmseqs2`.

    ```bash
    mmseqs2_msa --config_file <path_to_config>
    ```


**Prediction:**

- **predict_ensemble**: Run ensemble predictions.

    ```bash
    predict_ensemble --config_file <path_to_config>
    ```


- **fast_ensemble_init**: Initialize sample config file.

    ```bash
    fast_ensemble_init
    ```

**Analysis:**

- **rmsd_mode1d**: Analyze RMSD in 1D mode.

    ```bash
    rmsd_mode1d --config_file <path_to_config>
    ```

- **rmsd_mode2d**: Analyze RMSD in 2D mode.

    ```bash
    rmsd_mode2d --config_file <path_to_config>
    ```

- **tmscore_mode1d**: Analyze TM-score in 1D mode.

    ```bash
    tmscore_mode1d --config_file <path_to_config>
    ```

- **tmscore_mode2d**: Analyze TM-score in 2D mode.

    ```bash
    tmscore_mode2d --config_file <path_to_config>
    ```

- **pca_clustering**: Perform PCA clustering on the predicted structures.

    ```bash
    pca_clustering --config_file <path_to_config>
    ```

- **rmsf_plddt**: Calculate RMSF and pLDDT for the predicted structures.

    ```bash
    rmsf_plddt --config_file <path_to_config>
    ```

- **save_traj**: Save trajectories from the analysis.

    ```bash
    save_traj --config_file <path_to_config>
    ```
    

## Troubleshooting

- When running MSA generation or predictions, if you are getting a no such file or directory error, it is likely an issue with the AlphaFold2 Installation. Try reinstalling the package and ensure that the AF2 installation is on path.

- For issues with the qt platform, if the 'xcb' platform is found but cannot be initialized, try this command, or refer to [this](https://github.com/NVlabs/instant-ngp/discussions/300) github issue

```bash
sudo apt-get install libx11-xcb1 libxcb1 libxcb-glx0 \
libxcb-keysyms1 libxcb-image0 libxcb-shm0 libxcb-icccm4 \
libxcb-sync1 libxcb-xfixes0 libxcb-shape0 libxcb-randr0 \
libxcb-render-util0 libxcb-render0 libxcb-xinerama0 libxcb-xkb1 libxkbcommon-x11-0
```

## Download Sample Files

To get started quickly, download the sample files from the link below and add them to the root directory of the project:

[Download Sample Files](https://drive.google.com/drive/folders/1ev5HfWVyMTBw3FRtKWxYaswuaIXvC1FS?usp=drive_link)

## Contributing

We welcome contributions to FastEnsemble. If you would like to contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make your changes.
4. Commit your changes (`git commit -m 'Add some feature'`).
5. Push to the branch (`git push origin feature-branch`).
6. Create a new Pull Request.

Please ensure your code follows our coding standards and includes appropriate tests.

## License

FastEnsemble is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.


## Acknowledgements

Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1. 



