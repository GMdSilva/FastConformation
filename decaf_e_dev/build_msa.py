import warnings
from dataclasses import dataclass
from typing import Callable

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication, QLabel, QVBoxLayout, QWidget, QPushButton,
    QLineEdit, QCheckBox, QFileDialog, QComboBox
)
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'decaf_e_dev')))

from jackhmmer_msa import build_jackhmmer_msa
from mmseqs2_msa import run_mmseqs2_msa

class MSAOptionsWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        # MSA build method
        self.type_label = QLabel("MSA type:")
        self.type_dropdown = QComboBox()
        self.type_dropdown.addItems(["jackhmmer", "mmseqs2"])


        self.sequence_path_label = QLabel("Sequence Path:")
        self.sequence_path_input = QLineEdit()
        self.sequence_path_button = QPushButton("Select Sequence File")
        self.sequence_path_button.clicked.connect(self.select_sequence_path)
        
        self.tmp_dir_label = QLabel("Temporary Directory:")
        self.tmp_dir_input = QLineEdit("tmp")
        
        self.homooligomers_label = QLabel("Homooligomers:")
        self.homooligomers_input = QLineEdit("1")
        
        self.use_ramdisk_label = QLabel("Use RAM Disk:")
        self.use_ramdisk_checkbox = QCheckBox()
        self.use_ramdisk_checkbox.setChecked(True)
        
        layout.addWidget(self.type_label)
        layout.addWidget(self.type_dropdown)
        layout.addWidget(self.sequence_path_label)
        layout.addWidget(self.sequence_path_input)
        layout.addWidget(self.sequence_path_button)
        layout.addWidget(self.tmp_dir_label)
        layout.addWidget(self.tmp_dir_input)
        layout.addWidget(self.homooligomers_label)
        layout.addWidget(self.homooligomers_input)
        layout.addWidget(self.use_ramdisk_label)
        layout.addWidget(self.use_ramdisk_checkbox)
        
        self.setLayout(layout)
        self.setWindowTitle("MSA Options")

        # Run Button
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_msa)
        layout.addWidget(self.run_button)
        
    def select_sequence_path(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Sequence File", "", "FASTA Files (*.fasta);;All Files (*)", options=options)
        if file_path:
            self.sequence_path_input.setText(file_path)
    
    def run_msa(self):
        msa_type = self.type_dropdown.currentText()
        sequence_path = self.sequence_path_input.text()
        tmp_dir = self.tmp_dir_input.text()
        homooligomers = int(self.homooligomers_input.text())
        use_ramdisk = self.use_ramdisk_checkbox.isChecked()
        
        jobname = "msa_job"
        output_path = "msa_output"
        with open(sequence_path, 'r') as file:
            sequence_string = file.read()
        config = {
            "jobname": jobname,
            "sequence_path": sequence_path,
            "output_path": output_path,
            "use_ramdisk": use_ramdisk,
            "homooligomers": homooligomers,
            "tmp_dir": tmp_dir,
            "sequence_string": sequence_string
        }
        if msa_type == "jackhmmer":
            build_jackhmmer_msa(config)
        elif msa_type == "mmseqs2":
            run_mmseqs2_msa(config)

        print(f"MSA run with type: {msa_type}")