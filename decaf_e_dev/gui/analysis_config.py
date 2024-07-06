import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QStackedWidget, QCheckBox, QHBoxLayout, QSizePolicy
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon
from decaf_e_dev.gui.icons import Icons
from decaf_e_dev.rmsf_plddt import run_rmsf_analysis
from decaf_e_dev.rmsd_mode1d import run_rmsd_analysis
from decaf_e_dev.rmsd_mode2d import run_2d_rmsd_analysis
from decaf_e_dev.pca_clustering import run_pca_analysis
from decaf_e_dev.save_traj import run_trajectory_saving
from decaf_e_dev.tmscore_mode1d import run_tmscore_analysis
from decaf_e_dev.tmscore_mode2d import run_2d_tmscore_analysis
from decaf_e_dev.gui.widget_base import AnalysisWidgetBase, merge_configs

class AnalysisConfigWidget(QWidget):
    def __init__(self):
        super().__init__()
        main_layout = QVBoxLayout()
        
        # Top part for general analysis options and icon grid
        top_layout = QHBoxLayout()
        
        # General analysis options
        self.general_analysis_widget = GeneralAnalysisWidget()
        top_layout.addWidget(self.general_analysis_widget)
        
        # Icon grid
        self.icon_grid = Icons(self)
        self.icon_grid.addItems(ANALYSIS_CATEGORIES)
        self.icon_grid.itemClicked.connect(self._on_item_clicked)
        top_layout.addWidget(self.icon_grid)
        
        main_layout.addLayout(top_layout)
        
        # Stacked widget for analysis configuration panels
        self.analysis_stack = QStackedWidget()
        main_layout.addWidget(self.analysis_stack)
        
        # Plot display area
        self.plot_area = QLabel("Plot Area")
        self.plot_area.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.plot_area.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(self.plot_area)
        
        self.setLayout(main_layout)
        self.setWindowTitle("Analysis Configurations")

    def _on_item_clicked(self, item):
        name = item.text()
        getter = self.general_analysis_widget.get_general_options
        widget = ANALYSIS_CATEGORIES[name].widget(getter)
        
        self.analysis_stack.addWidget(widget)
        self.analysis_stack.setCurrentWidget(widget)
@dataclass
class AnalysisCategory:
    widget: Callable
    tool_tip: str = ""

ANALYSIS_CATEGORIES = {
    "RMSF Analysis": AnalysisCategory(
        widget=lambda getter: RMSFAnalysisWidget(getter),
        tool_tip="RMSF analysis configurations",
    ),
    "RMSD Analysis": AnalysisCategory(
        widget=lambda getter: RMSDAnalysisWidget(getter),
        tool_tip="RMSD analysis configurations",
    ),
    "RMSD_2D or TMSCORE_2D": AnalysisCategory(
        widget=lambda getter: RMSD2DWidget(getter),
        tool_tip="RMSD_2D or TMSCORE_2D configurations",
    ),
    "TMSCORE": AnalysisCategory(
        widget=lambda getter: TMSCOREWidget(getter),
        tool_tip="TMSCORE configurations",
    ),
    "PCA": AnalysisCategory(
        widget=lambda getter: PCAWidget(getter),
        tool_tip="PCA configurations",
    ),
    "Trajectory Saving": AnalysisCategory(
        widget=lambda getter: TrajectorySavingWidget(getter),
        tool_tip="Trajectory saving configurations",
    ),
}

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLineEdit, QPushButton, QFileDialog, QListWidget, QVBoxLayout
from PyQt5.QtCore import Qt
from dataclasses import dataclass
from typing import Callable

class GeneralAnalysisWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        
        self.jobname_input = QLineEdit("abl_wt")
        self.output_path_input = QLineEdit("/Users/fmgaleazzi/Downloads")
        self.output_path_button = QPushButton("Browse")
        self.predictions_path_input = QLineEdit("/Users/fmgaleazzi/decaf_e_dev/sample_predictions/abl_wt/predictions/alphafold2")
        self.predictions_path_button = QPushButton("Browse")
        self.engine_input = QLineEdit("alphafold2")
        self.align_range_input = QLineEdit("backbone")
        self.analysis_range_input = QLineEdit("backbone and resid 358-365")
        self.analysis_range_name_input = QLineEdit("aloop")
        self.ref1d_input = QLineEdit()
        self.starting_residue_input = QLineEdit("200")
        
        # Seq Pairs
        self.seq_pairs_layout = QVBoxLayout()
        self.add_seq_pair_button = QPushButton("Add Sequence Pair")
        self.add_seq_pair_button.clicked.connect(self.add_seq_pair)
        self.seq_pairs=[]
        layout.addRow("Job Name:", self.jobname_input)
        
        output_path_layout = QHBoxLayout()
        output_path_layout.addWidget(self.output_path_input)
        output_path_layout.addWidget(self.output_path_button)
        layout.addRow("Output Path:", output_path_layout)
        
        self.seq_pairs_layout.addWidget(self.add_seq_pair_button)
        layout.addRow("Sequence Pairs:", self.seq_pairs_layout)
        
        predictions_path_layout = QHBoxLayout()
        predictions_path_layout.addWidget(self.predictions_path_input)
        predictions_path_layout.addWidget(self.predictions_path_button)
        layout.addRow("Predictions Path:", predictions_path_layout)
        
        layout.addRow("Engine:", self.engine_input)
        layout.addRow("Align Range:", self.align_range_input)
        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        layout.addRow("Ref1D:", self.ref1d_input)
        layout.addRow("Starting Residue:", self.starting_residue_input)

        self.output_path_button.clicked.connect(self.select_output_path)
        self.predictions_path_button.clicked.connect(self.select_predictions_path)

        self.setLayout(layout)
        
    
    def select_output_path(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.output_path_input.setText(directory)

    def select_predictions_path(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.predictions_path_input.setText(directory)

    def get_general_options(self):
        return {
            "jobname": self.jobname_input.text(),
            "output_path": self.output_path_input.text(),
            "seq_pairs": [[int(seq_pair[0].text()), int(seq_pair[1].text())] for seq_pair in self.seq_pairs],
            "predictions_path": self.predictions_path_input.text(),
            "engine": self.engine_input.text(),
            "align_range": self.align_range_input.text(),
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "ref1d": self.ref1d_input.text(),
            "starting_residue": self.starting_residue_input.text(),
        }

    def add_seq_pair(self):
        seq_pair_layout = QHBoxLayout()
        seq_pair=[]
        seq1_input = QLineEdit()
        seq1_input.setPlaceholderText("Sequence 1")
        seq2_input = QLineEdit()
        seq2_input.setPlaceholderText("Sequence 2")
        seq_pair.append(seq1_input)
        seq_pair.append(seq2_input)
        self.seq_pairs.append(seq_pair)
        
        remove_button = QPushButton("Remove")
        remove_button.clicked.connect(lambda: self.remove_seq_pair(seq_pair_layout, seq_pair))

        seq_pair_layout.addWidget(seq1_input)
        seq_pair_layout.addWidget(seq2_input)
        seq_pair_layout.addWidget(remove_button)

        self.seq_pairs_layout.addLayout(seq_pair_layout)

    def remove_seq_pair(self, layout, seq_pair):
        # Properly remove and delete the layout and its widgets
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.seq_pairs_layout.removeItem(layout)
        self.seq_pairs.remove(seq_pair)
        layout.deleteLater()

class RMSFAnalysisWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.detect_mobile_checkbox = QCheckBox()
        self.detect_mobile_checkbox.setChecked(True)
        self.peak_width_input = QLineEdit("3")
        self.peak_prominence_input = QLineEdit("1")
        self.peak_height_input = QLineEdit("2")
        self.starting_residue_input = QLineEdit("200")
        
        layout.addRow("Detect Mobile:", self.detect_mobile_checkbox)
        layout.addRow("Peak Width:", self.peak_width_input)
        layout.addRow("Peak Prominence:", self.peak_prominence_input)
        layout.addRow("Peak Height:", self.peak_height_input)
        layout.addRow("Starting Residue:", self.starting_residue_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(lambda: self.run_analysis(general_options=True))
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.peak_width_input.text().isdigit():
            errors.append("Peak Width must be a number.")
        if not self.peak_prominence_input.text().isdigit():
            errors.append("Peak Prominence must be a number.")
        if not self.peak_height_input.text().isdigit():
            errors.append("Peak Height must be a number.")
        if not self.starting_residue_input.text().isdigit():
            errors.append("Starting Residue must be a number.")
        return errors

    def get_specific_options(self):
        return {
            "detect_mobile": self.detect_mobile_checkbox.isChecked(),
            "peak_width": int(self.peak_width_input.text()),
            "peak_prominence": int(self.peak_prominence_input.text()),
            "peak_height": int(self.peak_height_input.text()),
            "starting_residue": int(self.starting_residue_input.text())
        }

    def run_specific_analysis(self, config):
        run_rmsf_analysis(config)

class RMSDAnalysisWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.analysis_range_input = QLineEdit("backbone and resid 358-365")
        self.analysis_range_name_input = QLineEdit("aloop")
        
        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(lambda: self.run_analysis(general_options=True))
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.analysis_range_input.text():
            errors.append("Analysis Range cannot be empty.")
        if not self.analysis_range_name_input.text():
            errors.append("Analysis Range Name cannot be empty.")
        return errors

    def get_specific_options(self):
        return {
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text()
        }

    def run_specific_analysis(self, config):
        run_rmsd_analysis(config)
class RMSD2DWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        self.n_stdevs_input = QLineEdit("5")
        self.n_clusters_input = QLineEdit()

        layout.addRow("Mode Results:", self.mode_results_input)
        layout.addRow("Reference Structure 1 (Ref2D1):", self.ref2d1_input)
        layout.addRow("Reference Structure 2 (Ref2D2):", self.ref2d2_input)
        layout.addRow("Number of Standard Deviations (n_stdevs):", self.n_stdevs_input)
        layout.addRow("Number of Clusters (n_clusters):", self.n_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(lambda: self.run_analysis(general_options=True))
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.n_stdevs_input.text().isdigit():
            errors.append("Number of Standard Deviations (n_stdevs) must be a number.")
        return errors

    def get_specific_options(self):
        return {
            "mode_results": self.mode_results_input.text(),
            "ref2d1": self.ref2d1_input.text(),
            "ref2d2": self.ref2d2_input.text(),
            "n_stdevs": self.n_stdevs_input.text(),
            "n_clusters": self.n_clusters_input.text()
        }

    def run_specific_analysis(self, config):
        run_2d_tmscore_analysis(config)

class TMSCOREWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.slice_predictions_input = QLineEdit("backbone and resid 210-459")
        self.ref1_input = QLineEdit()

        layout.addRow("Slice Predictions:", self.slice_predictions_input)
        layout.addRow("Reference Structure (Ref1):", self.ref1_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.slice_predictions_input.text():
            errors.append("Slice Predictions cannot be empty.")
        return errors
    def get_specific_options(self):
        return {
            "slice_predictions": self.slice_predictions_input.text(),
            "ref1": self.ref1_input.text()
        }

    def run_specific_analysis(self, config):
        run_tmscore_analysis(config)
class TwoTMScoreWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        self.n_stdevs_input = QLineEdit("5")
        self.n_clusters_input = QLineEdit()

        layout.addRow("Mode Results:", self.mode_results_input)
        layout.addRow("Reference Structure 1 (Ref2D1):", self.ref2d1_input)
        layout.addRow("Reference Structure 2 (Ref2D2):", self.ref2d2_input)
        layout.addRow("Number of Standard Deviations (n_stdevs):", self.n_stdevs_input)
        layout.addRow("Number of Clusters (n_clusters):", self.n_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.mode_results_input.text():
            errors.append("Mode Results cannot be empty.")
        if not self.n_stdevs_input.text().isdigit():
            errors.append("Number of Standard Deviations (n_stdevs) must be a number.")
        return errors

    def get_specific_options(self):
        return {
            "mode_results": self.mode_results_input.text(),
            "ref2d1": self.ref2d1_input.text(),
            "ref2d2": self.ref2d2_input.text(),
            "n_stdevs": self.n_stdevs_input.text(),
            "n_clusters": self.n_clusters_input.text()
        }

    def run_specific_analysis(self, config):
        run_2d_tmscore_analysis(config)

class PCAWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.align_range_input = QLineEdit()
        self.analysis_range_input = QLineEdit()
        self.analysis_range_name_input = QLineEdit()
        self.n_pca_clusters_input = QLineEdit("3")

        layout.addRow("Align Range:", self.align_range_input)
        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        layout.addRow("Number of PCA Clusters:", self.n_pca_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.n_pca_clusters_input.text().isdigit():
            errors.append("Number of PCA Clusters must be a number.")
        return errors

    def get_specific_options(self):
        return {
            "align_range": self.align_range_input.text(),
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "n_pca_clusters": self.n_pca_clusters_input.text()
        }

    def run_specific_analysis(self, config):
        run_pca_analysis(config)

class TrajectorySavingWidget(AnalysisWidgetBase):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.analysis_range_input = QLineEdit()
        self.analysis_range_name_input = QLineEdit()
        self.reorder_input = QLineEdit("rmsd_1d")
        self.traj_format_input = QLineEdit("pdb")

        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        layout.addRow("Reorder By:", self.reorder_input)
        layout.addRow("Trajectory Format:", self.traj_format_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        errors = []
        if not self.analysis_range_input.text():
            errors.append("Analysis Range cannot be empty.")
        if not self.analysis_range_name_input.text():
            errors.append("Analysis Range Name cannot be empty.")
        if not self.reorder_input.text():
            errors.append("Reorder By cannot be empty.")
        if not self.traj_format_input.text():
            errors.append("Trajectory Format cannot be empty.")
        return errors

    def get_specific_options(self):
        return {
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "reorder": self.reorder_input.text(),
            "traj_format": self.traj_format_input.text()
        }

    def run_specific_analysis(self, config):
        run_trajectory_saving(config)