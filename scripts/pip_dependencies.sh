#!/bin/bash

conda install -c biocore hmmer

echo "Installing Python dependencies..."
pip install -e .

echo "Upgrading jax and jaxlib..."
pip install --upgrade jax jaxlib