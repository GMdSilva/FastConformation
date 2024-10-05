import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
import os
from PyQt5.QtWidgets import (
    QApplication, QScrollArea, QVBoxLayout, QWidget, QLabel, QStackedWidget, QCheckBox, QHBoxLayout, QMessageBox
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon
from fast_conformation.gui.icons import Icons
from fast_conformation.rmsf_plddt import run_rmsf_analysis
from fast_conformation.rmsd_mode1d import run_rmsd_analysis
from fast_conformation.rmsd_mode2d import run_2d_rmsd_analysis
from fast_conformation.pca_clustering import run_pca_analysis
from fast_conformation.save_traj import run_trajectory_saving
from fast_conformation.tmscore_mode1d import run_tmscore_analysis
from fast_conformation.tmscore_mode2d import run_2d_tmscore_analysis
from fast_conformation.gui.widget_base import AnalysisWidgetBase, merge_configs
from fast_conformation.gui.plot_widget import PlotWidget
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLineEdit, QPushButton, QFileDialog, QListWidget, QVBoxLayout
from PyQt5.QtCore import Qt
from dataclasses import dataclass
from typing import Callable


class AnalysisConfigWidget(QWidget):
    """
    The AnalysisConfigWidget class is responsible for providing a user interface 
    to configure various analysis options. It contains general analysis options 
    and specific configurations for different types of analysis (RMSF, RMSD, TMScore, etc.).

    Attributes:
        job_manager: A manager object that handles job submissions and execution.
    """
    def __init__(self, job_manager):
        """
        Initialize the AnalysisConfigWidget with a job manager.

        Args:
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__()
        self.job_manager = job_manager

        # Main layout
        main_layout = QVBoxLayout()
        
        # Top part for general analysis options and icon grid
        top_layout = QHBoxLayout()
        self.setStyleSheet("""
            QToolBar {
                background-color: #333333;
                color: white;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton {
                background-color: #555555;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 0 5px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #666666;
            }
            QPushButton:pressed {
                background-color: #777777;
            }
        """)
        
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
        
        # Create a scroll area and set the main layout as its widget
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_content.setLayout(main_layout)
        scroll_area.setWidget(scroll_content)
        
        # Set the scroll area as the central widget
        layout = QVBoxLayout()
        layout.addWidget(scroll_area)
        self.setLayout(layout)
        self.setWindowTitle("Analysis Configurations")

    def _on_item_clicked(self, item):
        """
        Handle the event when an item is clicked in the icon grid. It loads the 
        corresponding analysis widget in the stacked widget.

        Args:
            item: The clicked item in the icon grid.
        """
        name = item.text()
        getter = self.general_analysis_widget.get_general_options
        widget = ANALYSIS_CATEGORIES[name].widget(getter, self.job_manager)
        widget.setStyleSheet("""
            QToolBar {
                background-color: #333333;
                color: white;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton {
                background-color: #555555;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 0 5px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #666666;
            }
            QPushButton:pressed {
                background-color: #777777;
            }
        """)
        self.analysis_stack.addWidget(widget)
        self.analysis_stack.setCurrentWidget(widget)

    def show_plot(self, name, plot_widget):
        """
        Display the plot for the given analysis.

        Args:
            name: The name of the analysis.
            plot_widget: The widget containing the plot to be displayed.
        """
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget()
        print(f'parent_widget {parent}')
        parent.show_plot(name, plot_widget)

@dataclass
class AnalysisCategory:
    """
    Data class to represent an analysis category.

    Attributes:
        widget: A callable that returns the widget associated with the analysis.
        tool_tip: A string that provides a description of the analysis category.
    """
    widget: Callable
    tool_tip: str = ""

ANALYSIS_CATEGORIES = {
    "RMSF Analysis": AnalysisCategory(
        widget=lambda getter, job_manager: RMSFAnalysisWidget(getter, job_manager),
        tool_tip="RMSF analysis configurations",
    ),
    "RMSD Analysis": AnalysisCategory(
        widget=lambda getter, job_manager: RMSDAnalysisWidget(getter, job_manager),
        tool_tip="RMSD analysis configurations",
    ),
    "RMSD_2D": AnalysisCategory(
        widget=lambda getter, job_manager: RMSD2DWidget(getter, job_manager),
        tool_tip="RMSD_2D configurations",
    ),
    "TMSCORE": AnalysisCategory(
        widget=lambda getter, job_manager: TMSCOREWidget(getter, job_manager),
        tool_tip="TMSCORE configurations",
    ),
    "TMSCORE_2D": AnalysisCategory(
        widget=lambda getter, job_manager: TwoTMScoreWidget(getter, job_manager),
        tool_tip="TMSCORE_2D configurations",
    ),
    "PCA": AnalysisCategory(
        widget=lambda getter, job_manager: PCAWidget(getter, job_manager),
        tool_tip="PCA configurations",
    ),
    "Trajectory Saving": AnalysisCategory(
        widget=lambda getter, job_manager: TrajectorySavingWidget(getter, job_manager),
        tool_tip="Trajectory saving configurations",
    ),
}

class GeneralAnalysisWidget(QWidget):
    """
    GeneralAnalysisWidget provides the user interface for setting general analysis 
    options like job name, output path, sequence pairs, and other configurations.

    Methods:
        clear_seq_pairs: Clears the list of sequence pairs.
        select_output_path: Opens a file dialog to select the output path directory.
        select_predictions_path: Opens a file dialog to select the predictions path directory.
        auto_detect_sequence_pairs: Automatically detects sequence pairs from the selected directory.
        get_general_options: Returns a dictionary of the current general analysis options.
        add_seq_pair: Adds a new sequence pair to the list.
        remove_seq_pair: Removes a sequence pair from the list.
    """
    def __init__(self):
        """
        Initialize the GeneralAnalysisWidget with input fields and layout.
        """
        super().__init__()
        layout = QFormLayout()
        
        self.jobname_input = QLineEdit("abl_wt")
        self.output_path_input = QLineEdit("")
        self.output_path_button = QPushButton("Browse")
        self.predictions_path_input = QLineEdit("")
        self.predictions_path_button = QPushButton("Browse")
        self.engine_input = QLineEdit("alphafold2")
        self.align_range_input = QLineEdit("backbone")
        self.analysis_range_input = QLineEdit("backbone and resid 358-365")
        self.analysis_range_name_input = QLineEdit("aloop")
        self.ref1d_input = QLineEdit()
        self.starting_residue_input = QLineEdit("200")
        
        # Seq Pairs
        self.seq_pairs_layout = QVBoxLayout()
        self.add_seq_pair_button = QPushButton("Add Pair:")
        self.add_seq_pair_button.clicked.connect(lambda: self.add_seq_pair("", ""))
        self.seq_pairs = []
        
        layout.addRow("Job Name:", self.jobname_input)
        
        output_path_layout = QHBoxLayout()
        output_path_layout.addWidget(self.output_path_input)
        output_path_layout.addWidget(self.output_path_button)
        layout.addRow("Output Path:", output_path_layout)
        
        instructions_label = QLabel("Select the path to AF2 predictions which contains subfolders for each subsampling condition.")
        instructions_label.setWordWrap(True)  # Ensure the text wraps if it's too long
        layout.addRow(instructions_label)

        predictions_path_layout = QHBoxLayout()
        predictions_path_layout.addWidget(self.predictions_path_input)
        predictions_path_layout.addWidget(self.predictions_path_button)
        layout.addRow("Predictions Path:", predictions_path_layout)
        
        self.seq_pairs_layout.addWidget(self.add_seq_pair_button)
        layout.addRow("max_seq:extra_seq pairs:", self.seq_pairs_layout)
        
        layout.addRow("Engine:", self.engine_input)
        layout.addRow("Align Range:", self.align_range_input)
        layout.addRow("Analysis Range:", self.analysis_range_input)
        layout.addRow("Analysis Range Name:", self.analysis_range_name_input)
        layout.addRow("Ref1D:", self.ref1d_input)
        layout.addRow("Starting Residue:", self.starting_residue_input)

        self.output_path_button.clicked.connect(self.select_output_path)
        self.predictions_path_button.clicked.connect(self.select_predictions_path)

        self.setLayout(layout)
        
    def clear_seq_pairs(self):
        """
        Clears all sequence pairs from the layout.
        """
        for i in range(len(self.seq_pairs) - 1, -1, -1):
            layout, seq_pair = self.seq_pairs.pop(i)
            while layout.count():
                child = layout.takeAt(0)
                if child.widget():
                    child.widget().deleteLater()
            self.seq_pairs_layout.removeItem(layout)
            layout.deleteLater()

    def select_output_path(self):
        """
        Opens a file dialog to allow the user to select an output directory.
        """
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.output_path_input.setText(directory)

    def select_predictions_path(self):
        """
        Opens a file dialog to allow the user to select a directory containing AF2 predictions.
        """
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        if directory:
            self.predictions_path_input.setText(directory)
            self.auto_detect_sequence_pairs(directory)

    def auto_detect_sequence_pairs(self, predictions_path):
        """
        Automatically detects sequence pairs based on the structure of the predictions path directory.

        Args:
            predictions_path: The path to the directory containing predictions.
        """
        try:
            self.clear_seq_pairs()
            for subdir in os.listdir(predictions_path):
                subdir_path = os.path.join(predictions_path, subdir)
                if os.path.isdir(subdir_path):
                    parts = subdir.split('_')
                    if len(parts) >= 3 and parts[-2].isdigit() and parts[-1].isdigit():
                        max_seq = parts[-2]
                        extra_seq = parts[-1]
                        self.add_seq_pair(max_seq, extra_seq)
                    else:
                        current_pair = subdir.split('_')[-2:]
                        if all(part.isdigit() for part in current_pair):
                            max_seq, extra_seq = current_pair
                            self.add_seq_pair(max_seq, extra_seq)
                            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to auto-detect sequence pairs: {e}")
    
    def get_general_options(self):
        """
        Retrieves the current general options set in the widget.

        Returns:
            A dictionary containing the general options.
        """
        return {
            "jobname": self.jobname_input.text(),
            "output_path": self.output_path_input.text(),
            "seq_pairs": [[int(seq_pair[0].text()), int(seq_pair[1].text())] for layout, seq_pair in self.seq_pairs],
            "predictions_path": self.predictions_path_input.text(),
            "engine": self.engine_input.text(),
            "align_range": self.align_range_input.text(),
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "ref1d": self.ref1d_input.text(),
            "starting_residue": self.starting_residue_input.text(),
        }

    def add_seq_pair(self, seq1='', seq2=''):
        """
        Adds a new sequence pair to the widget.

        Args:
            seq1: The first sequence in the pair.
            seq2: The second sequence in the pair.
        """
        seq_pair_layout = QHBoxLayout()
        seq_pair = []
        seq1_input = QLineEdit(seq1)
        seq1_input.setPlaceholderText("Sequence 1")
        seq2_input = QLineEdit(seq2)
        seq2_input.setPlaceholderText("Sequence 2")
        seq_pair.append(seq1_input)
        seq_pair.append(seq2_input)
        
        remove_button = QPushButton("Remove")
        remove_button.clicked.connect(lambda: self.remove_seq_pair(seq_pair_layout, seq_pair))

        seq_pair_layout.addWidget(seq1_input)
        seq_pair_layout.addWidget(seq2_input)
        seq_pair_layout.addWidget(remove_button)

        self.seq_pairs_layout.insertLayout(self.seq_pairs_layout.count() - 1, seq_pair_layout)
        self.seq_pairs.append((seq_pair_layout, seq_pair))

    def remove_seq_pair(self, layout, seq_pair):
        """
        Removes a sequence pair from the widget.

        Args:
            layout: The layout containing the sequence pair widgets.
            seq_pair: The sequence pair to be removed.
        """
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.seq_pairs_layout.removeItem(layout)
        self.seq_pairs.remove((layout, seq_pair))
        layout.deleteLater()

class RMSFAnalysisWidget(AnalysisWidgetBase):
    """
    Widget to configure and run RMSF (Root Mean Square Fluctuation) analysis.

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for RMSF analysis.
        run_specific_analysis: Runs the RMSF analysis with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the RMSFAnalysisWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
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

        self.layout_rmsf=layout

        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        """
        Validates the inputs provided by the user for RMSF analysis.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
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
        """
        Retrieves the specific options for RMSF analysis.

        Returns:
            A dictionary containing the specific options for RMSF analysis.
        """
        return {
            "detect_mobile": self.detect_mobile_checkbox.isChecked(),
            "peak_width": int(self.peak_width_input.text()),
            "peak_prominence": int(self.peak_prominence_input.text()),
            "peak_height": int(self.peak_height_input.text()),
            "starting_residue": int(self.starting_residue_input.text())
        }

    def run_specific_analysis(self, config):
        """
        Runs the RMSF analysis with the given configuration.

        Args:
            config: The configuration options for the RMSF analysis.
        """
        self.plot_widget = PlotWidget(self)
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
        run_rmsf_analysis(config, self.plot_widget)
        parent.show_plot("RMSF Analysis", self.plot_widget)

class RMSDAnalysisWidget(AnalysisWidgetBase):
    """
    Widget to configure and run RMSD (Root Mean Square Deviation) analysis.

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for RMSD analysis.
        run_specific_analysis: Runs the RMSD analysis with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the RMSDAnalysisWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
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

    def validate_inputs(self):
        """
        Validates the inputs provided by the user for RMSD analysis.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
        errors = []
        if not self.analysis_range_input.text():
            errors.append("Analysis Range cannot be empty.")
        if not self.analysis_range_name_input.text():
            errors.append("Analysis Range Name cannot be empty.")
        return errors

    def get_specific_options(self):
        """
        Retrieves the specific options for RMSD analysis.

        Returns:
            A dictionary containing the specific options for RMSD analysis.
        """
        return {
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text()
        }

    def run_specific_analysis(self, config):
        """
        Runs the RMSD analysis with the given configuration.

        Args:
            config: The configuration options for the RMSD analysis.
        """
        self.plot_widget = PlotWidget(self)
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
        run_rmsd_analysis(config, self.plot_widget)
        parent.show_plot("RMSD Analysis", self.plot_widget)

class RMSD2DWidget(AnalysisWidgetBase):
    """
    Widget to configure and run 2D RMSD (Root Mean Square Deviation) analysis.

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for 2D RMSD analysis.
        run_specific_analysis: Runs the 2D RMSD analysis with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the RMSD2DWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        self.n_stdevs_input = QLineEdit("5")
        self.n_clusters_input = QLineEdit()

        layout.addRow("RMSD1D analysis results csv filepath:", self.mode_results_input)
        layout.addRow("Reference Structure 1 (Ref2D1):", self.ref2d1_input)
        layout.addRow("Reference Structure 2 (Ref2D2):", self.ref2d2_input)
        layout.addRow("Number of Standard Deviations (n_stdevs):", self.n_stdevs_input)
        layout.addRow("Number of Clusters (n_clusters):", self.n_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        """
        Validates the inputs provided by the user for 2D RMSD analysis.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
        errors = []
        if not self.n_stdevs_input.text().isdigit():
            errors.append("Number of Standard Deviations (n_stdevs) must be a number.")
        return errors

    def get_specific_options(self):
        """
        Retrieves the specific options for 2D RMSD analysis.

        Returns:
            A dictionary containing the specific options for 2D RMSD analysis.
        """
        return {
            "mode_results": self.mode_results_input.text(),
            "ref2d1": self.ref2d1_input.text(),
            "ref2d2": self.ref2d2_input.text(),
            "n_stdevs": self.n_stdevs_input.text(),
            "n_clusters": self.n_clusters_input.text()
        }

    def run_specific_analysis(self, config):
        """
        Runs the 2D RMSD analysis with the given configuration.

        Args:
            config: The configuration options for the 2D RMSD analysis.
        """
        self.plot_widget = PlotWidget(self)
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
        run_2d_rmsd_analysis(config, self.plot_widget)
        parent.show_plot("TMScore-2D Analysis", self.plot_widget)

class TMSCOREWidget(AnalysisWidgetBase):
    """
    Widget to configure and run TMScore analysis.

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for TMScore analysis.
        run_specific_analysis: Runs the TMScore analysis with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the TMSCOREWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
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
        """
        Validates the inputs provided by the user for TMScore analysis.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
        errors = []
        return errors
    
    def get_specific_options(self):
        """
        Retrieves the specific options for TMScore analysis.

        Returns:
            A dictionary containing the specific options for TMScore analysis.
        """
        return {
            "slice_predictions": self.slice_predictions_input.text(),
            "ref1": self.ref1_input.text()
        }

    def run_specific_analysis(self, config):
        """
        Runs the TMScore analysis with the given configuration.

        Args:
            config: The configuration options for the TMScore analysis.
        """
        self.plot_widget = PlotWidget(self)
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
        run_tmscore_analysis(config, self.plot_widget)
        parent.show_plot("TMScore Analysis", self.plot_widget)

class TwoTMScoreWidget(AnalysisWidgetBase):
    """
    Widget to configure and run 2D TMScore analysis.

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for 2D TMScore analysis.
        run_specific_analysis: Runs the 2D TMScore analysis with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the TwoTMScoreWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.mode_results_input = QLineEdit()
        self.ref2d1_input = QLineEdit()
        self.ref2d2_input = QLineEdit()
        self.n_stdevs_input = QLineEdit("5")
        self.n_clusters_input = QLineEdit()

        layout.addRow("TMScore1D analysis results csv filepath:", self.mode_results_input)
        layout.addRow("Reference Structure 1 (Ref2D1):", self.ref2d1_input)
        layout.addRow("Reference Structure 2 (Ref2D2):", self.ref2d2_input)
        layout.addRow("Number of Standard Deviations (n_stdevs):", self.n_stdevs_input)
        layout.addRow("Number of Clusters (n_clusters):", self.n_clusters_input)
        
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_button)
        
        self.setLayout(layout)

    def validate_inputs(self):
        """
        Validates the inputs provided by the user for 2D TMScore analysis.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
        errors = []
        if not self.mode_results_input.text():
            errors.append("Mode Results cannot be empty.")
        if not self.n_stdevs_input.text().isdigit():
            errors.append("Number of Standard Deviations (n_stdevs) must be a number.")
        return errors

    def get_specific_options(self):
        """
        Retrieves the specific options for 2D TMScore analysis.

        Returns:
            A dictionary containing the specific options for 2D TMScore analysis.
        """
        return {
            "mode_results": self.mode_results_input.text(),
            "ref2d1": self.ref2d1_input.text(),
            "ref2d2": self.ref2d2_input.text(),
            "n_stdevs": self.n_stdevs_input.text(),
            "n_clusters": self.n_clusters_input.text()
        }

    def run_specific_analysis(self, config):
        """
        Runs the 2D TMScore analysis with the given configuration.

        Args:
            config: The configuration options for the 2D TMScore analysis.
        """
        self.plot_widget = PlotWidget(self)
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
        run_2d_tmscore_analysis(config, self.plot_widget)
        parent.show_plot("TMScore-2D Analysis", self.plot_widget)

class PCAWidget(AnalysisWidgetBase):
    """
    Widget to configure and run PCA (Principal Component Analysis).

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for PCA analysis.
        run_specific_analysis: Runs the PCA analysis with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the PCAWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
        self.general_options_getter = general_options_getter
        layout = QFormLayout()
        
        self.align_range_input = QLineEdit("backbone")
        self.analysis_range_input = QLineEdit("backbone and resid 358-365")
        self.analysis_range_name_input = QLineEdit("aloop")
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
        """
        Validates the inputs provided by the user for PCA analysis.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
        errors = []
        if not self.align_range_input.text():
            errors.append("Align Range cannot be empty.")
        if not self.analysis_range_input.text():
            errors.append("Analysis Range cannot be empty.")
        if not self.analysis_range_name_input.text():
            errors.append("Analysis Range Name cannot be empty.")
        if not self.n_pca_clusters_input.text().isdigit():
            errors.append("Number of PCA Clusters must be a number.")
        return errors
    
    def get_specific_options(self):
        """
        Retrieves the specific options for PCA analysis.

        Returns:
            A dictionary containing the specific options for PCA analysis.
        """
        return {
            "align_range": self.align_range_input.text(),
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "n_pca_clusters": self.n_pca_clusters_input.text()
        }

    def run_specific_analysis(self, config):
        """
        Runs the PCA analysis with the given configuration.

        Args:
            config: The configuration options for the PCA analysis.
        """
        self.plot_widget = PlotWidget(self)
        parent=self.parentWidget().parentWidget().parentWidget().parentWidget().parentWidget()
        run_pca_analysis(config, self.plot_widget)
        parent.show_plot("TMScore-2D Analysis", self.plot_widget)

class TrajectorySavingWidget(AnalysisWidgetBase):
    """
    Widget to configure and run trajectory saving.

    Methods:
        validate_inputs: Validates the inputs provided by the user.
        get_specific_options: Retrieves the specific options for trajectory saving.
        run_specific_analysis: Runs the trajectory saving process with the given configuration.
    """
    def __init__(self, general_options_getter, job_manager):
        """
        Initialize the TrajectorySavingWidget with general options getter and job manager.

        Args:
            general_options_getter: Callable to retrieve general analysis options.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(job_manager)
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
        """
        Validates the inputs provided by the user for trajectory saving.

        Returns:
            A list of errors if any input is invalid, otherwise an empty list.
        """
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
        """
        Retrieves the specific options for trajectory saving.

        Returns:
            A dictionary containing the specific options for trajectory saving.
        """
        return {
            "analysis_range": self.analysis_range_input.text(),
            "analysis_range_name": self.analysis_range_name_input.text(),
            "reorder": self.reorder_input.text(),
            "traj_format": self.traj_format_input.text()
        }

    def run_specific_analysis(self, config):
        """
        Runs the trajectory saving process with the given configuration.

        Args:
            config: The configuration options for trajectory saving.
        """
        run_trajectory_saving(config)