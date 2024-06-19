# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Common utilities for data pipeline tools."""
import contextlib
import shutil
import tempfile
import time
from typing import Optional
import subprocess
from getpass import getpass
import os

from absl import logging


@contextlib.contextmanager
def tmpdir_manager(base_dir: Optional[str] = None):
  """Context manager that deletes a temporary directory on exit."""
  tmpdir = tempfile.mkdtemp(dir=base_dir)
  try:
    yield tmpdir
  finally:
    shutil.rmtree(tmpdir, ignore_errors=True)


@contextlib.contextmanager
def timing(msg: str):
  logging.info('Started %s', msg)
  tic = time.time()
  yield
  toc = time.time()
  logging.info('Finished %s in %.3f seconds', msg, toc - tic)


def create_ram_disk():
  password = getpass("Enter your sudo password to create Ram Disk for running HMMER faster: ")

  # Create the ramdisk directory
  command_mkdir = "sudo mkdir -m 777 --parents /tmp/ramdisk"
  subprocess.run(['sudo', '-S'] + command_mkdir.split(), input=password.encode())

  # Mount the ramdisk
  command_mount = "sudo mount -t tmpfs -o size=9G ramdisk /tmp/ramdisk"
  subprocess.run(['sudo', '-S'] + command_mount.split(), input=password.encode())


def read_fasta(file_path):
  if not os.path.exists(file_path):
    raise FileNotFoundError(f"FASTA file {file_path} not found.")

  sequences = {}
  with open(file_path, 'r') as file:
    sequence_name = ''
    sequence_data = ''
    valid_fasta = False
    for line in file:
      line = line.strip()
      if line.startswith('>'):
        if sequence_name:
          if not sequence_data:
            raise ValueError(f"Sequence data for {sequence_name} is missing.")
          sequences[sequence_name] = sequence_data
          sequence_data = ''
        sequence_name = line[1:]  # Remove the '>' character
        valid_fasta = True
      elif line:
        if not sequence_name:
          raise ValueError("FASTA file is missing a sequence header before sequence data.")
        sequence_data += line
    if sequence_name:  # Add the last sequence to the dictionary
      if not sequence_data:
        raise ValueError(f"Sequence data for {sequence_name} is missing.")
      sequences[sequence_name] = sequence_data

  if not valid_fasta:
    raise ValueError("No valid FASTA format detected in the file.")

  return sequences


def save_dict_to_fasta(seq_dict, output_path, jobname):
  for seq_name, seq in seq_dict.items():
    with open(f'{output_path}/{jobname}/target_seq/{jobname}.fasta', 'w') as file:
      file.write(f">{jobname}\n")
      file.write(f"{seq}\n")
      return ## TODO this is a hack to get just the first sequence, replace with something more elegant


def create_directory(path):
  if os.path.exists(path):
    shutil.rmtree(path)
  os.makedirs(path)
  print(f"Directory '{path}' created successfully.")