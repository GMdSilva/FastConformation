#!/bin/bash

# Set the download URL for the TM-align binary
TALIGN_URL="https://zhanggroup.org/TM-align/TMalign.tar.gz"

# Set the target installation directory (default to /usr/local/bin)
INSTALL_DIR="/usr/local/bin"

# Download the TM-align binary
echo "Downloading TM-align from $TALIGN_URL..."
wget -O TMalign.tar.gz $TALIGN_URL

# Check if download was successful
if [ $? -ne 0 ]; then
  echo "Error: Failed to download TM-align. Please check your internet connection and the URL."
  exit 1
fi

# Extract the downloaded file
echo "Extracting TM-align..."
tar -xzvf TMalign.tar.gz

# Check if extraction was successful
if [ $? -ne 0 ]; then
  echo "Error: Failed to extract TM-align."
  exit 1
fi

# Move the TMalign executable to the installation directory
echo "Installing TM-align to $INSTALL_DIR..."
sudo mv TMalign $INSTALL_DIR

# Check if the move was successful
if [ $? -ne 0 ]; then
  echo "Error: Failed to move TM-align to $INSTALL_DIR. You may need to run this script as root or with sudo."
  exit 1
fi

# Make sure the executable has the correct permissions
echo "Setting executable permissions for TM-align..."
sudo chmod +x $INSTALL_DIR/TMalign

# Verify the installation
echo "Verifying the installation..."
if command -v TMalign &> /dev/null; then
  echo "TM-align has been successfully installed and is available in your PATH."
else
  echo "Error: TM-align is not found in your PATH. Please check the installation directory."
  exit 1
fi

# Clean up the downloaded file
rm TMalign.tar.gz

echo "Installation completed successfully!"
