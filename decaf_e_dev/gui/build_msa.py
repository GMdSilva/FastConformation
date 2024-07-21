from decaf_e_dev.jackhmmer_msa import build_jackhmmer_msa
from decaf_e_dev.mmseqs2_msa import run_mmseqs2_msa
from decaf_e_dev.gui.widget_base import AnalysisWidgetBase, merge_configs
from PyQt5.QtWidgets import QFileDialog, QLabel, QVBoxLayout, QComboBox, QLineEdit, QPushButton, QCheckBox, QMessageBox
from PyQt5.QtCore import Qt

class MSAOptionsWidget(AnalysisWidgetBase):
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
        self.use_ramdisk_checkbox.setChecked(False)

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
        self.run_button.clicked.connect(lambda: self.run_analysis(general_options=False))
        layout.addWidget(self.run_button)
        
    def select_sequence_path(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Sequence File", "", "FASTA Files (*.fasta);;All Files (*)", options=options)
        if file_path:
            self.sequence_path_input.setText(file_path)

    def validate_inputs(self):
        errors = []
        if not self.sequence_path_input.text():
            errors.append("Sequence Path cannot be empty.")
        if not self.tmp_dir_input.text():
            errors.append("Temporary Directory cannot be empty.")
        if not self.homooligomers_input.text().isdigit():
            errors.append("Homooligomers must be a number.")
        return errors

    def get_specific_options(self):
        sequence_path = self.sequence_path_input.text()
        with open(sequence_path, 'r') as file:
            sequence_string = file.read()
        
        return {
            "jobname": "msa_job",
            "sequence_path": sequence_path,
            "output_path": "msa_output",
            "use_ramdisk": self.use_ramdisk_checkbox.isChecked(),
            "homooligomers": int(self.homooligomers_input.text()),
            "tmp_dir": self.tmp_dir_input.text(),
            "sequence_string": sequence_string,
            "msa_type": self.type_dropdown.currentText()
        }
    

    def run_specific_analysis(self, config):
        msa_type = config["msa_type"]
        if msa_type == "jackhmmer":
            build_jackhmmer_msa(config)
        elif msa_type == "mmseqs2":
            run_mmseqs2_msa(config)

        print(f"MSA run with type: {msa_type}")
