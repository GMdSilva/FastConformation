#!/bin/bash

# Define the Miniconda installer URL
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
MINICONDA_INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
MINICONDA_DIR="$HOME/miniconda3"

# Function to prompt user for confirmation
confirm() {
    read -p "$1 (y/n): " choice
    case "$choice" in
        y|Y ) return 0;;
        n|N ) return 1;;
        * ) echo "Invalid input"; confirm "$1";;
    esac
}

# Check if Miniconda is already installed
if [ -d "$MINICONDA_DIR" ]; then
    echo "Miniconda is already installed."
else
    # Download and install Miniconda
    if confirm "Do you want to download and install Miniconda?"; then
        if [ -f "$MINICONDA_INSTALLER" ]; then
            echo "Miniconda installer already downloaded."
        else
            echo "Downloading Miniconda installer..."
            wget $MINICONDA_URL -O $MINICONDA_INSTALLER
        fi

        echo "Installing Miniconda..."
        bash $MINICONDA_INSTALLER -b
    else
        echo "Skipping Miniconda installation."
    fi
fi

# Check if Conda is initialized
if conda -V &> /dev/null; then
    echo "Conda is already initialized."
else
    if confirm "Do you want to initialize Conda?"; then
        echo "Initializing Conda..."
        $MINICONDA_DIR/bin/conda init
        echo "Refreshing shell..."
        source ~/.bashrc
    else
        echo "Skipping Conda initialization."
    fi
fi

# Check if Mamba is installed
if conda list | grep mamba &> /dev/null; then
    echo "Mamba is already installed."
else
    if confirm "Do you want to install Mamba?"; then
        echo "Installing Mamba..."
        $MINICONDA_DIR/bin/conda install mamba -n base -c conda-forge -y
    else
        echo "Skipping Mamba installation."
    fi
fi

# Check if the environment already exists
if conda env list | grep decaf_e_dev &> /dev/null; then
    echo "Environment 'decaf_e_dev' already exists."
else
    if confirm "Do you want to create the Mamba environment 'decaf_e_dev' with Python 3.9?"; then
        echo "Creating Mamba environment 'decaf_e_dev' with Python 3.9..."
        $MINICONDA_DIR/bin/mamba create -n decaf_e_dev python=3.9 -y
    else
        echo "Skipping environment creation."
    fi
fi

# Check if the environment is already active
if [ "$CONDA_DEFAULT_ENV" = "decaf_e_dev" ]; then
    echo "Environment 'decaf_e_dev' is already active."
else
    if confirm "Do you want to activate the 'decaf_e_dev' environment?"; then
        echo "Activating the 'decaf_e_dev' environment..."
        source $MINICONDA_DIR/bin/activate decaf_e_dev
        echo "Environment 'decaf_e_dev' is now active."
    else
        echo "Skipping environment activation."
    fi
fi

# Install dependencies
if confirm "Do you want to install dependencies?"; then
    echo "Installing TM-SCore."
    bash scripts/install_tmscore.sh

    echo "Installing local-colabfold."
    bash scripts/install_colabfold.sh

    echo "Installing hmmer."
    bash scripts/install_hmmer_aptget.sh

    echo "Installing Python dependencies..."
    pip install -r requirements.txt
    pip install -e .

    echo "Upgrading jax and jaxlib..."
    pip install --upgrade jax jaxlib
else
    echo "Skipping dependencies installation."
fi
