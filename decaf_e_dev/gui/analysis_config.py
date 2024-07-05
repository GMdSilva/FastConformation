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
        
        self.jobname_input = QLineEdit()
        self.output_path_input = QLineEdit()
        self.output_path_button = QPushButton("Browse")
        self.seq_pairs_input = QLineEdit()
        self.predictions_path_input = QLineEdit()
        self.predictions_path_button = QPushButton("Browse")
        self.engine_input = QLineEdit()
        self.align_range_input = QLineEdit("backbone")
        self.analysis_range_input = QLineEdit()
        self.analysis_range_name_input = QLineEdit()
        self.ref1d_input = QLineEdit()
        self.starting_residue_input = QLineEdit("200")
        
        layout.addRow("Job Name:", self.jobname_input)
        
        output_path_layout = QHBoxLayout()
        output_path_layout.addWidget(self.output_path_input)
        output_path_layout.addWidget(self.output_path_button)
        layout.addRow("Output Path:", output_path_layout)
        
        layout.addRow("Sequence Pairs:", self.seq_pairs_input)
        
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
            "seq_pairs": self.seq_pairs_input.text(),
            "predictions_path": self.predictions_path_input.text(),
            "engine": self.engine_input.text(),
            "align_range": self.align_range_input.text(),
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "ref1d": self.ref1d_input.text(),
            "starting_residue": self.starting_residue_input.text(),
        }

class RMSFAnalysisWidget(QWidget):
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
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def run_analysis(self):
        general_options = self.general_options_getter()

        detect_mobile = self.detect_mobile_checkbox.isChecked()
        peak_width = int(self.peak_width_input.text())
        peak_prominence = int(self.peak_prominence_input.text())
        peak_height = int(self.peak_height_input.text())
        
        self.plot_window = PlotWindow()
        # Run the analysis using the general and specific options
        run_rmsf_analysis(
            config_file=None,
            jobname=general_options["jobname"],
            output_path=general_options["output_path"],
            seq_pairs=general_options["seq_pairs"],
            predictions_path=general_options["predictions_path"],
            engine=general_options["engine"],
            align_range=general_options["align_range"],
            detect_mobile=detect_mobile,
            peak_width=peak_width,
            peak_prominence=peak_prominence,
            peak_height=peak_height,
            starting_residue=general_options["starting_residue"]
        )

class RMSDAnalysisWidget(QWidget):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.analysis_range_input = QLineEdit("backbone and resid 358-365")
        self.analysis_range_name_input = QLineEdit("aloop")
        
        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        
        self.setLayout(layout)

    def run_analysis(self):
        general_options=self.general_options_getter()
        
        # Get specific options for RMSD analysis
        ref1d = self.ref1d_input.text()
        analysis_range = self.analysis_range_input.text()
        analysis_range_name = self.analysis_range_name_input.text()

        # Run the analysis using the general and specific options
        run_rmsd_analysis(
            config_file=None,
            jobname=general_options["jobname"],
            output_path=general_options["output_path"],
            seq_pairs=general_options["seq_pairs"],
            predictions_path=general_options["predictions_path"],
            engine=general_options["engine"],
            align_range=general_options["align_range"],
            analysis_range=analysis_range,
            analysis_range_name=analysis_range_name,
            ref1d=ref1d,
            starting_residue=general_options["starting_residue"]
        )

class RMSD2DWidget(QWidget):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        self.n_stdevs_input = QLineEdit()
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

    def run_analysis(self):
        general_options = self.general_options_getter()

        # Get specific options for 2D RMSD analysis
        mode_results = self.mode_results_input.text()
        ref2d1 = self.ref2d1_input.text()
        ref2d2 = self.ref2d2_input.text()
        n_stdevs = self.n_stdevs_input.text()
        n_clusters = self.n_clusters_input.text()

        # Run the analysis using the general and specific options
        run_2d_rmsd_analysis(
            config_file=None,
            jobname=general_options["jobname"],
            output_path=general_options["output_path"],
            mode_results=mode_results,
            seq_pairs=general_options["seq_pairs"],
            predictions_path=general_options["predictions_path"],
            engine=general_options["engine"],
            align_range=general_options["align_range"],
            analysis_range=general_options["analysis_range"],
            analysis_range_name=general_options["analysis_range_name"],
            ref2d1=ref2d1,
            ref2d2=ref2d2,
            n_stdevs=n_stdevs,
            n_clusters=n_clusters,
            starting_residue=general_options["starting_residue"]
        )
class TMSCOREWidget(QWidget):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.slice_predictions_input = QLineEdit()
        self.ref1_input = QLineEdit()

        layout.addRow("Slice Predictions:", self.slice_predictions_input)
        layout.addRow("Reference Structure (Ref1):", self.ref1_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def run_analysis(self):
        general_options = self.general_options_getter()

        # Get specific options for TM-Score analysis
        slice_predictions = self.slice_predictions_input.text()
        ref1 = self.ref1_input.text()

        # Run the analysis using the general and specific options
        run_tmscore_analysis(
            config_file=None,
            output_path=general_options["output_path"],
            predictions_path=general_options["predictions_path"],
            jobname=general_options["jobname"],
            seq_pairs=general_options["seq_pairs"],
            starting_residue=general_options["starting_residue"],
            slice_predictions=slice_predictions,
            ref1=ref1,
            engine=general_options["engine"]
        )

class TwoTMScoreWidget(QWidget):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        self.n_stdevs_input = QLineEdit()
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

    def run_analysis(self):
        general_options = self.general_options_getter()

        # Get specific options for 2D TM-Score analysis
        mode_results = self.mode_results_input.text()
        ref2d1 = self.ref2d1_input.text()
        ref2d2 = self.ref2d2_input.text()
        n_stdevs = self.n_stdevs_input.text()
        n_clusters = self.n_clusters_input.text()

        # Run the analysis using the general and specific options
        run_2d_tmscore_analysis(
            config_file=None,
            output_path=general_options["output_path"],
            predictions_path=general_options["predictions_path"],
            mode_results=mode_results,
            jobname=general_options["jobname"],
            seq_pairs=general_options["seq_pairs"],
            starting_residue=general_options["starting_residue"],
            slice_predictions=general_options["slice_predictions"],
            ref2d1=ref2d1,
            ref2d2=ref2d2,
            engine=general_options["engine"],
            n_stdevs=n_stdevs,
            n_clusters=n_clusters
        )

class PCAWidget(QWidget):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.align_range_input = QLineEdit()
        self.analysis_range_input = QLineEdit()
        self.analysis_range_name_input = QLineEdit()
        self.n_pca_clusters_input = QLineEdit()

        layout.addRow("Align Range:", self.align_range_input)
        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        layout.addRow("Number of PCA Clusters:", self.n_pca_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def run_analysis(self):
        general_options = self.general_options_getter()

        # Get specific options for PCA analysis
        align_range = self.align_range_input.text()
        analysis_range = self.analysis_range_input.text()
        analysis_range_name = self.analysis_range_name_input.text()
        n_pca_clusters = self.n_pca_clusters_input.text()

        # Run the analysis using the general and specific options
        run_pca_analysis(
            config_file=None,
            predictions_path=general_options["predictions_path"],
            output_path=general_options["output_path"],
            seq_pairs=general_options["seq_pairs"],
            jobname=general_options["jobname"],
            align_range=align_range,
            analysis_range=analysis_range,
            analysis_range_name=analysis_range_name,
            engine=general_options["engine"],
            n_pca_clusters=n_pca_clusters,
            starting_residue=general_options["starting_residue"]
        )

class TrajectorySavingWidget(QWidget):
    def __init__(self, general_options_getter):
        super().__init__()
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.analysis_range_input = QLineEdit()
        self.analysis_range_name_input = QLineEdit()
        self.reorder_input = QLineEdit()
        self.traj_format_input = QLineEdit()

        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        layout.addRow("Reorder By:", self.reorder_input)
        layout.addRow("Trajectory Format:", self.traj_format_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def run_analysis(self):
        general_options = self.general_options_getter()

        # Get specific options for trajectory saving
        analysis_range = self.analysis_range_input.text()
        analysis_range_name = self.analysis_range_name_input.text()
        reorder = self.reorder_input.text()
        traj_format = self.traj_format_input.text()

        # Run the analysis using the general and specific options
        run_trajectory_saving(
            config_file=None,
            output_path=general_options["output_path"],
            predictions_path=general_options["predictions_path"],
            jobname=general_options["jobname"],
            seq_pairs=general_options["seq_pairs"],
            starting_residue=general_options["starting_residue"],
            analysis_range=analysis_range,
            analysis_range_name=analysis_range_name,
            reorder=reorder,
            traj_format=traj_format,
            engine=general_options["engine"]
        )
