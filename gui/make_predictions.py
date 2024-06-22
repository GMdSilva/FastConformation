import warnings
from dataclasses import dataclass
from typing import Callable

from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QApplication, QLabel, QVBoxLayout, QWidget, QPushButton,
    QLineEdit, QCheckBox, QFileDialog, QComboBox, QListWidget, QGridLayout, QHBoxLayout
)

class MakePredictionsWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_ui()

    def init_ui(self):
        layout = QGridLayout()

        # Engine
        self.engine_label = QLabel("Engine:")
        self.engine_dropdown = QComboBox()
        self.engine_dropdown.addItems(["alphafold2", "other_engine1", "other_engine2"])

        # MSA Path
        self.msa_path_label = QLabel("MSA Path:")
        self.msa_path_input = QLineEdit()
        self.msa_path_button = QPushButton("Select MSA File")
        self.msa_path_button.clicked.connect(self.select_msa_path)

        # MSA From
        self.msa_from_label = QLabel("MSA From:")
        self.msa_from_dropdown = QComboBox()
        self.msa_from_dropdown.addItems(["mmseqs2", "other_source1", "other_source2"])

        # Seq Pairs
        self.seq_pairs_label = QLabel("Sequence Pairs:")
        self.seq_pairs_layout = QVBoxLayout()
        self.add_seq_pair_button = QPushButton("Add Sequence Pair")
        self.add_seq_pair_button.clicked.connect(self.add_seq_pair)

        # Seeds
        self.seeds_label = QLabel("Seeds:")
        self.seeds_input = QLineEdit("10")

        # Platform
        self.platform_label = QLabel("Platform:")
        self.platform_dropdown = QComboBox()
        self.platform_dropdown.addItems(["cpu", "gpu"])

        # Save All
        self.save_all_label = QLabel("Save All:")
        self.save_all_checkbox = QCheckBox()
        self.save_all_checkbox.setChecked(False)

        # Models
        self.models_label = QLabel("Models:")
        self.models_input = QLineEdit("[1, 2, 3, 4, 5]")

        # Recycles
        self.recycles_label = QLabel("Recycles:")
        self.recycles_input = QLineEdit("4")

        # Subset MSA To
        self.subset_msa_to_label = QLabel("Subset MSA To:")
        self.subset_msa_to_input = QLineEdit("")

        # Adding widgets to the layout
        layout.addWidget(self.engine_label, 0, 0)
        layout.addWidget(self.engine_dropdown, 0, 1)
        layout.addWidget(self.msa_path_label, 1, 0)
        layout.addWidget(self.msa_path_input, 1, 1)
        layout.addWidget(self.msa_path_button, 1, 2)
        layout.addWidget(self.msa_from_label, 2, 0)
        layout.addWidget(self.msa_from_dropdown, 2, 1)
        layout.addWidget(self.seq_pairs_label, 3, 0)
        layout.addLayout(self.seq_pairs_layout, 3, 1, 1, 2)
        layout.addWidget(self.add_seq_pair_button, 4, 1)
        layout.addWidget(self.seeds_label, 5, 0)
        layout.addWidget(self.seeds_input, 5, 1)
        layout.addWidget(self.platform_label, 6, 0)
        layout.addWidget(self.platform_dropdown, 6, 1)
        layout.addWidget(self.save_all_label, 7, 0)
        layout.addWidget(self.save_all_checkbox, 7, 1)
        layout.addWidget(self.models_label, 8, 0)
        layout.addWidget(self.models_input, 8, 1)
        layout.addWidget(self.recycles_label, 9, 0)
        layout.addWidget(self.recycles_input, 9, 1)
        layout.addWidget(self.subset_msa_to_label, 10, 0)
        layout.addWidget(self.subset_msa_to_input, 10, 1)

        self.setLayout(layout)
        self.setWindowTitle("Advanced MSA Options")

    def select_msa_path(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Select MSA File", "", "MSA Files (*.m3a);;All Files (*)", options=options)
        if file_path:
            self.msa_path_input.setText(file_path)

    def add_seq_pair(self):
        seq_pair_layout = QHBoxLayout()

        seq1_input = QLineEdit()
        seq1_input.setPlaceholderText("Sequence 1")
        seq2_input = QLineEdit()
        seq2_input.setPlaceholderText("Sequence 2")

        remove_button = QPushButton("Remove")
        remove_button.clicked.connect(lambda: self.remove_seq_pair(seq_pair_layout))

        seq_pair_layout.addWidget(seq1_input)
        seq_pair_layout.addWidget(seq2_input)
        seq_pair_layout.addWidget(remove_button)

        self.seq_pairs_layout.addLayout(seq_pair_layout)

    def remove_seq_pair(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
