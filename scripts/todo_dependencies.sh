
# Update the package lists
sudo apt-get update

# Install g++
sudo apt-get install -y g++

# Download TMscore.cpp
sudo apt install tm-align

#!/bin/bash

conda install -c biocore hmmer

echo "Installing Python dependencies..."
pip install -e .

echo "Upgrading jax and jaxlib..."
pip install --upgrade jax jaxlib

#!/bin/bash

# Step 1: Download the installation script
sudo apt install --quiet --yes hmmer

echo "Installation of hmmer complete."