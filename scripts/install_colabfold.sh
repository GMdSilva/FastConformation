#!/bin/bash

# Step 1: Download the installation script
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh

# Step 2: Execute the installation script
bash install_colabbatch_linux.sh

# Step 3: Add localcolabfold to PATH
INSTALL_DIR="$(pwd)/localcolabfold/colabfold-conda/bin"
PROFILE_FILE="$HOME/.bashrc"

# Check if the PATH is already present in the profile file
if ! grep -q "$INSTALL_DIR" "$PROFILE_FILE"; then
    echo "Adding $INSTALL_DIR to PATH in $PROFILE_FILE"
    echo "export PATH=\"$INSTALL_DIR:\$PATH\"" >> "$PROFILE_FILE"
    source "$PROFILE_FILE"
else
    echo "$INSTALL_DIR is already in PATH"
fi

echo "Installation of colabfold complete."
