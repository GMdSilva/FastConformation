import warnings
from dataclasses import dataclass
from typing import Callable

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication, QLabel, QVBoxLayout, QWidget, QPushButton,
    QLineEdit, QCheckBox, QFileDialog
)

class MSAOptionsWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

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

    def select_sequence_path(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Sequence File", "", "FASTA Files (*.fasta);;All Files (*)", options=options)
        if file_path:
            self.sequence_path_input.setText(file_path)
