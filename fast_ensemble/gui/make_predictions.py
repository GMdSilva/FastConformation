from PyQt5.QtWidgets import QApplication, QWidget, QGridLayout, QLabel, QComboBox, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog, QCheckBox
import sys
from fast_ensemble.predict_ensemble import run_ensemble_prediction
from fast_ensemble.predict_ensemble import run_ensemble_prediction
from fast_ensemble.gui.widget_base import AnalysisWidgetBase, merge_configs
from PyQt5.QtWidgets import QVBoxLayout, QGridLayout, QLabel, QComboBox, QLineEdit, QPushButton, QCheckBox, QFileDialog, QHBoxLayout

class MakePredictionsWidget(AnalysisWidgetBase):
    """
    The MakePredictionsWidget class provides a user interface for configuring 
    and running predictions using different MSA (Multiple Sequence Alignment) 
    options and settings.

    Methods:
        init_ui: Initializes the user interface components.
        select_msa_path: Opens a file dialog to select the MSA file.
        select_output_path: Opens a file dialog to select the output directory.
        add_seq_pair: Adds a new sequence pair input to the interface.
        remove_seq_pair: Removes an existing sequence pair input.
        validate_inputs: Validates the user inputs to ensure they are correct.
        get_specific_options: Retrieves the specific options set by the user.
        run_analysis: Validates inputs, merges configurations, and starts the prediction job.
        get_seq_pairs: Retrieves the sequence pairs input by the user.
    """
    def __init__(self, job_manager, general_options_getter=None, *args, **kwargs):
        """
        Initialize the MakePredictionsWidget with a job manager and optional general options getter.

        Args:
            job_manager: The manager responsible for handling job execution.
            general_options_getter: An optional callable to retrieve general analysis options.
            *args: Additional arguments to pass to the parent class.
            **kwargs: Additional keyword arguments to pass to the parent class.
        """
        super().__init__(job_manager, *args, **kwargs)
        self.general_options_getter = general_options_getter
        self.init_ui()

    def init_ui(self):
        """
        Initializes the user interface components and layout for the widget.
        """
        layout = QGridLayout()

        # Engine
        self.engine_label = QLabel("Engine:")
        self.engine_dropdown = QComboBox()
        self.engine_dropdown.addItems(["alphafold2"])
        self.setStyleSheet("""
            QToolBar {
                background-color: #333333;
                color: white;
                padding: 10px;
            }
            QPushButton {
                background-color: #555555;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 0 5px;
            }
            QPushButton:hover {
                background-color: #666666;
            }
            QPushButton:pressed {
                background-color: #777777;
            }
        """)
        
        self.job_name_label = QLabel("MSA Path:")
        self.job_name_input = QLineEdit()

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

        # Output Path
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
        layout.addWidget(self.job_name_label, 12, 0)
        layout.addWidget(self.job_name_input, 12, 1)        

        self.setLayout(layout)
        self.setWindowTitle("Advanced MSA Options")
        self.add_seq_pair(seq1="256", seq2="512")

        # Run Button
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(lambda: self.run_analysis())
        self.output_path_button.clicked.connect(self.select_output_path)
        layout.addWidget(self.run_button, 12, 1)
        
    def select_msa_path(self):
        """
        Opens a file dialog to allow the user to select the MSA file.
        """
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(self, "Select MSA File", "", "MSA Files (*.a3m);;All Files (*)", options=options)
        if file_path:
            self.msa_path_input.setText(file_path)
    
    def select_output_path(self):
        """
        Opens a file dialog to allow the user to select the output directory.
        """
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.output_path_input.setText(directory)

    def add_seq_pair(self, seq1="", seq2=""):
        """
        Adds a new sequence pair input to the interface.

        Args:
            seq1: The first sequence in the pair.
            seq2: The second sequence in the pair.
        """
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
        """
        Removes an existing sequence pair input from the interface.

        Args:
            layout: The layout containing the sequence pair to be removed.
        """
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.seq_pairs_layout.removeItem(layout)
        layout.deleteLater()

    def validate_inputs(self):
        """
        Validates the user inputs to ensure they are correct.

        Returns:
            A list of error messages if validation fails, otherwise an empty list.
        """
        errors = []
        if not self.msa_path_input.text():
            errors.append("MSA Path cannot be empty.")
        if not self.seeds_input.text().isdigit():
            errors.append("Seeds must be a number.")
        if self.subset_msa_to_input.text() and not self.subset_msa_to_input.text().isdigit():
            errors.append("Subset MSA To must be a number.")
        return errors

    def get_specific_options(self):
        """
        Retrieves the specific options set by the user.

        Returns:
            A dictionary containing the options for the prediction job.
        """
        return {
            'msa_path': self.msa_path_input.text(),
            'output_path': self.output_path_input.text(),
            'jobname': self.job_name_input.text(),
            'seq_pairs': self.get_seq_pairs(),
            'seeds': int(self.seeds_input.text()),
            'save_all': self.save_all_checkbox.isChecked(),
            'platform': self.platform_dropdown.currentText(),
            'subset_msa_to': int(self.subset_msa_to_input.text()) if self.subset_msa_to_input.text() else None,
            'msa_from': self.msa_from_dropdown.currentText()
        }

    def run_analysis(self):
        """
        Validates inputs, merges configurations, and starts the prediction job.
        """
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
        self.show_info_message(f"Job {config['jobname']} started.")

    def get_seq_pairs(self):
        """
        Retrieves the sequence pairs input by the user.

        Returns:
            A list of sequence pairs, each represented as a list of two integers.
        """
        seq_pairs = []
        for i in range(self.seq_pairs_layout.count()):
            layout = self.seq_pairs_layout.itemAt(i).layout()
            seq1 = layout.itemAt(0).widget().text()
            seq2 = layout.itemAt(1).widget().text()
            seq_pairs.append([int(seq1), int(seq2)])
        return seq_pairs
