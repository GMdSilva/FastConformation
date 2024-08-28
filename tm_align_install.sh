#!/bin/bash

# Script to install TM-score on a Linux system

# Define TM-score URL and installation directory
TMSCORE_URL="https://zhanggroup.org/TM-score/TMscore.cpp"
INSTALL_DIR="$HOME/tmscore"

# Create the installation directory
mkdir -p "$INSTALL_DIR"

# Download TM-score source code
echo "Downloading TM-score source code..."
curl -o "$INSTALL_DIR/TMscore.cpp" "$TMSCORE_URL"

# Navigate to the installation directory
cd "$INSTALL_DIR"

# Compile TM-score
echo "Compiling TM-score..."
g++ -O3 -ffast-math -lm -o TMscore TMscore.cpp

# Check if compilation was successful
if [ -f "TMscore" ]; then
    echo "TM-score successfully installed!"

    # Add TM-score to PATH
    echo "Adding TM-score to your PATH..."
    echo "export PATH=\"\$PATH:$INSTALL_DIR\"" >> ~/.bashrc
    source ~/.bashrc

    echo "TM-score installation complete! You can use TMscore command from anywhere in the terminal."
else
    echo "Error: TM-score compilation failed. Please check the output for details."
fi
