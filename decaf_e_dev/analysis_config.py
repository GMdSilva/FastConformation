import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QStackedWidget, QListWidget, QListWidgetItem, QHBoxLayout, QSizePolicy
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon
from icons import Icons

class AnalysisConfigWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)

        general_analysis_widget = GeneralAnalysisWidget()
        layout.addWidget(general_analysis_widget)

        self.icon_grid = Icons(self)
        self.icon_grid.addItems(ANALYSIS_CATEGORIES)
        self.icon_grid.itemClicked.connect(self._on_item_clicked)
        layout.addWidget(self.icon_grid)

        self.setLayout(layout)
        self.setWindowTitle("Analysis Configurations")

    def _on_item_clicked(self, item):
        name = item.text()
        widget = ANALYSIS_CATEGORIES[name].widget()

        # Clear any existing widgets before adding new one
        while self.layout().count() > 3:
            item = self.layout().takeAt(3)
            item.widget().deleteLater()

        self.layout().addWidget(widget)
@dataclass
class AnalysisCategory:
    widget: Callable
    tool_tip: str = ""

ANALYSIS_CATEGORIES = {
    "RMSF Analysis": AnalysisCategory(
        widget=lambda: RMSFAnalysisWidget(),
        tool_tip="RMSF analysis configurations",
    ),
    "RMSD Analysis": AnalysisCategory(
        widget=lambda: RMSDAnalysisWidget(),
        tool_tip="RMSD analysis configurations",
    ),
    "RMSD_2D or TMSCORE_2D": AnalysisCategory(
        widget=lambda: RMSD2DOrTMSCORE2DWidget(),
        tool_tip="RMSD_2D or TMSCORE_2D configurations",
    ),
    "TMSCORE": AnalysisCategory(
        widget=lambda: TMSCOREWidget(),
        tool_tip="TMSCORE configurations",
    ),
    "PCA": AnalysisCategory(
        widget=lambda: PCAWidget(),
        tool_tip="PCA configurations",
    ),
    "Trajectory Saving": AnalysisCategory(
        widget=lambda: TrajectorySavingWidget(),
        tool_tip="Trajectory saving configurations",
    ),
}
from qtpy.QtWidgets import QPushButton, QTextEdit, QLineEdit, QCheckBox, QFormLayout

class GeneralAnalysisWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        
        self.align_range_input = QLineEdit("backbone")
        self.ref1d_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        
        layout.addRow("Align Range:", self.align_range_input)
        layout.addRow("Ref1D:", self.ref1d_input)
        layout.addRow("Ref2D1:", self.ref2d1_input)
        layout.addRow("Ref2D2:", self.ref2d2_input)

        self.setLayout(layout)

class RMSFAnalysisWidget(QWidget):
    def __init__(self):
        super().__init__()
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
        # Simulate analysis and display results
        results = f"Detect Mobile: {self.detect_mobile_checkbox.isChecked()}\n"
        results += f"Peak Width: {self.peak_width_input.text()}\n"
        results += f"Peak Prominence: {self.peak_prominence_input.text()}\n"
        results += f"Peak Height: {self.peak_height_input.text()}\n"
        results += f"Starting Residue: {self.starting_residue_input.text()}\n"


class RMSDAnalysisWidget(QWidget):
    def __init__(self):
        super().__init__()
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
        # Simulate analysis and display results
        results = f"Analysis Range: {self.analysis_range_input.text()}\n"
        results += f"Analysis Range Name: {self.analysis_range_name_input.text()}\n"

class RMSD2DOrTMSCORE2DWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.n_stdevs_input = QLineEdit("5")
        
        layout.addRow("Mode Results:", self.mode_results_input)
        layout.addRow("N Stdevs:", self.n_stdevs_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def run_analysis(self):
        # Simulate analysis and display results
        results = f"Mode Results: {self.mode_results_input.text()}\n"
        results += f"N Stdevs: {self.n_stdevs_input.text()}\n"


class TMSCOREWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        
        self.slice_predictions_input = QLineEdit("backbone and resid 210-459")
        
        layout.addRow("Slice Predictions:", self.slice_predictions_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)

        
        self.setLayout(layout)

    def run_analysis(self):
        # Simulate analysis and display results
        results = f"Slice Predictions: {self.slice_predictions_input.text()}\n"


class PCAWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        
        self.n_pca_clusters_input = QLineEdit("3")
        
        layout.addRow("N PCA Clusters:", self.n_pca_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)

        
        self.setLayout(layout)

    def run_analysis(self):
        # Simulate analysis and display results
        results = f"N PCA Clusters: {self.n_pca_clusters_input.text()}\n"


class TrajectorySavingWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        
        self.reorder_input = QLineEdit("rmsd_1d")
        self.traj_format_input = QLineEdit("pdb")
        
        layout.addRow("Reorder:", self.reorder_input)
        layout.addRow("Trajectory Format:", self.traj_format_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)

        
        self.setLayout(layout)

    def run_analysis(self):
        # Simulate analysis and display results
        results = f"Reorder: {self.reorder_input.text()}\n"
        results += f"Trajectory Format: {self.traj_format_input.text()}\n"
 
