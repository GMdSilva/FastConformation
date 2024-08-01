from PyQt5.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QComboBox, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog, QCheckBox
import sys
from decaf_e_dev.predict_ensemble import run_ensemble_prediction
from decaf_e_dev.predict_ensemble import run_ensemble_prediction
from decaf_e_dev.gui.widget_base import AnalysisWidgetBase, merge_configs

class MakePredictionsWidget(AnalysisWidgetBase):
    def __init__(self, job_manager, general_options_getter=None, *args, **kwargs):
        super().__init__(job_manager, *args, **kwargs)
        self.general_options_getter = general_options_getter
        self.init_ui()

    def init_ui(self):
        layout = QGridLayout()

        # Engine
        self.engine_label = QLabel("Engine:")
        self.engine_dropdown = QComboBox()
        self.engine_dropdown.addItems(["alphafold2"])

        # MSA Path
        self.msa_path_label = QLabel("MSA Path:")
        self.msa_path_input = QLineEdit()
        self.msa_path_button = QPushButton("Select MSA File")
        self.msa_path_button.clicked.connect(self.select_msa_path)

        # MSA From
        self.msa_from_label = QLabel("MSA From:")
        self.msa_from_dropdown = QComboBox()
        self.msa_from_dropdown.addItems(["mmseqs2", "jackhmmer"])

        # Seq Pairs
        self.seq_pairs_label = QLabel("max_seq:extra_seq pairs:")
        self.seq_pairs_layout = QVBoxLayout()
        self.add_seq_pair_button = QPushButton("Add Pair")
        self.add_seq_pair_button.clicked.connect(lambda: self.add_seq_pair(seq1="", seq2=""))

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

        # Subset MSA To
        self.subset_msa_to_label = QLabel("Max MSA Depth:")
        self.subset_msa_to_input = QLineEdit("")
        self.output_path_label = QLabel("Output Path:")
        self.output_path_input = QLineEdit("")
        self.output_path_button = QPushButton("Browse")
        
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
        layout.addWidget(self.subset_msa_to_label, 10, 0)
        layout.addWidget(self.subset_msa_to_input, 10, 1)
        layout.addWidget(self.output_path_label, 11, 0)
        layout.addWidget(self.output_path_input, 11, 1)
        layout.addWidget(self.output_path_button, 11, 2)

        self.setLayout(layout)
        self.setWindowTitle("Advanced MSA Options")
        self.add_seq_pair(seq1="256", seq2="512")
        # Run Button
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(lambda: self.run_analysis())
        self.output_path_button.clicked.connect(self.select_output_path)
        layout.addWidget(self.run_button, 12, 1)
        
    def select_msa_path(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Select MSA File", "", "MSA Files (*.a3m);;All Files (*)", options=options)
        if file_path:
            self.msa_path_input.setText(file_path)
    
    def select_output_path(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.output_path_input.setText(directory)

    def add_seq_pair(self, seq1="", seq2=""):
        seq_pair_layout = QHBoxLayout()

        seq1_input = QLineEdit(seq1)
        seq1_input.setPlaceholderText("Sequence 1")
        seq2_input = QLineEdit(seq2)
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
        self.seq_pairs_layout.removeItem(layout)
        layout.deleteLater()

    def validate_inputs(self):
        errors = []
        if not self.msa_path_input.text():
            errors.append("MSA Path cannot be empty.")
        if not self.seeds_input.text().isdigit():
            errors.append("Seeds must be a number.")
        if self.subset_msa_to_input.text() and not self.subset_msa_to_input.text().isdigit():
            errors.append("Subset MSA To must be a number.")
        return errors

    def get_specific_options(self):
        return {
            'msa_path': self.msa_path_input.text(),
            'output_path': self.output_path_input.text(),  # Set this as needed
            'jobname': 'jobname',  # Set this as needed
            'seq_pairs': self.get_seq_pairs(),
            'seeds': int(self.seeds_input.text()),
            'save_all': self.save_all_checkbox.isChecked(),
            'platform': self.platform_dropdown.currentText(),
            'subset_msa_to': int(self.subset_msa_to_input.text()) if self.subset_msa_to_input.text() else None,
            'msa_from': self.msa_from_dropdown.currentText()
        }

    def run_analysis(self):
        errors = self.validate_inputs()
        if errors:
            self.show_error_message(errors)
            return

        try:
            g_options = self.general_options_getter()
        except Exception:
            g_options = {}

        specific_options = self.get_specific_options()
        config = merge_configs(g_options, specific_options)

        job_id = self.job_manager.run_job(run_ensemble_prediction, (config,), config['jobname'])
        self.show_info_message(f"Job {job_id} started.")

    def get_seq_pairs(self):
        seq_pairs = []
        for i in range(self.seq_pairs_layout.count()):
            layout = self.seq_pairs_layout.itemAt(i).layout()
            seq1 = layout.itemAt(0).widget().text()
            seq2 = layout.itemAt(1).widget().text()
            seq_pairs.append([int(seq1), int(seq2)])
        return seq_pairs
