#!/bin/bash

# Path to the README.md file
README_PATH="README.md"

# Path to the requirements.txt file
REQ_PATH="requirements.txt"

# Destination directory
DEST_DIR="docs"

# Destination path for the README.md file in the docs directory
DEST_README_PATH="$DEST_DIR/README.md"

# Destination path for the requirements.txt file in the docs directory
DEST_REQ_PATH="$DEST_DIR/requirements.txt"

# Check if README.md exists at the specified path
if [ ! -f "$README_PATH" ]; then
  echo "README.md not found at $README_PATH"
  exit 1
fi

# Check if requirements.txt exists at the specified path
if [ ! -f "$REQ_PATH" ]; then
  echo "requirements.txt not found at $REQ_PATH"
  exit 1
fi

# Copy the README.md into the docs directory, overwriting if it exists
cp "$README_PATH" "$DEST_README_PATH"

# Copy the requirements.txt into the docs directory, overwriting if it exists
cp "$REQ_PATH" "$DEST_REQ_PATH"

# Navigate to the docs directory
cd "$DEST_DIR" || exit

# Clean previous builds
make clean

# Build the HTML documentation
make html

# Notify completion
echo "Documentation build complete."
