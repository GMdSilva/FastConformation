#!/bin/bash

# Define the Miniconda installer URL
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

# Function to prompt user for confirmation
confirm() {
    read -p "$1 (y/n): " choice
    case "$choice" in
        y|Y ) return 0;;
        n|N ) return 1;;
        * ) echo "Invalid input"; confirm "$1";;
    esac
}

# Download and install Miniconda
if confirm "Do you want to download and install Miniconda?"; then
    echo "Downloading Miniconda installer..."
    wget $MINICONDA_URL -O Miniconda3-latest-Linux-x86_64.sh

    echo "Installing Miniconda..."
    bash Miniconda3-latest-Linux-x86_64.sh -b
else
    echo "Skipping Miniconda installation."
fi

# Initialize Conda
if confirm "Do you want to initialize Conda?"; then
    echo "Initializing Conda..."
    ~/miniconda3/bin/conda init
    echo "Refreshing shell..."
    source ~/.bashrc
else
    echo "Skipping Conda initialization."
fi

# Install Mamba
if confirm "Do you want to install Mamba?"; then
    echo "Installing Mamba..."
    ~/miniconda3/bin/conda install mamba -n base -c conda-forge -y
else
    echo "Skipping Mamba installation."
fi

# Create Mamba environment with Python 3.9
if confirm "Do you want to create the Mamba environment 'decaf_e_dev' with Python 3.9?"; then
    echo "Creating Mamba environment 'decaf_e_dev' with Python 3.9..."
    ~/miniconda3/bin/mamba create -n decaf_e_dev python=3.9 -y
else
    echo "Skipping environment creation."
fi

# Activate the environment
if confirm "Do you want to activate the 'decaf_e_dev' environment?"; then
    echo "Activating the 'decaf_e_dev' environment..."
    source ~/miniconda3/bin/activate decaf_e_dev
    echo "Environment 'decaf_e_dev' is now active."
else
    echo "Skipping environment activation."
fi

# Install dependencies
if confirm "Do you want to install dependencies?"; then
    echo "Installing TM-SCore."
    bash scripts/install_tmscore.sh

    echo "Installing local-colabfold."
    bash scripts/install_colabfold.sh

    echo "Installing hmmer."
    bash scripts.install_hmmer_aptget.sh

    pip install -r requirements.txt
    pip install -e .

    pip install --upgrade jax jaxlib
else
    echo "Skipping dependencies installation."
fi
