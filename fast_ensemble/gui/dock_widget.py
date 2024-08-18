from dataclasses import dataclass
from typing import Callable
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QLabel, QVBoxLayout, QWidget, QDockWidget, QPushButton, QHBoxLayout, QSizePolicy
)
from qtpy.QtGui import QFont
from fast_ensemble.gui.icons import Icons
from fast_ensemble.gui.build_msa import MSAOptionsWidget
from fast_ensemble.gui.make_predictions import MakePredictionsWidget
from fast_ensemble.gui.analysis_config import AnalysisConfigWidget
from fast_ensemble.gui.job_manager import JobStatusPage, JobManager
from dataclasses import dataclass
from typing import Callable
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QSizePolicy, QDockWidget, Qt

@dataclass
class Category:
    """
    A data class that represents a category in the main widget.

    Attributes:
        widget: A callable that returns a QWidget, taking a JobManager as an argument.
        tool_tip: A string representing the tooltip text for the category.
    """
    widget: Callable[[JobManager], QWidget]
    tool_tip: str = ""

CATEGORIES = {
    "Build MSA": Category(
        widget=lambda job_manager: MSAOptionsWidget(job_manager),
        tool_tip="Select parameters to build MSA",
    ),
    "Make Predictions": Category(
        widget=lambda job_manager: MakePredictionsWidget(job_manager),
        tool_tip="Select parameters to make predictions",
    ),
    "Analysis": Category(
        widget=lambda job_manager: AnalysisConfigWidget(job_manager),
        tool_tip="Select parameters to analyze results",
    ),
}


class MainWidget(QWidget):
    """
    MainWidget serves as the central widget for managing different tasks such as 
    building MSA, making predictions, and analyzing results. It provides an 
    interface for selecting tasks and displaying the appropriate configuration widgets.

    Methods:
        create_dock_widget: Creates a new dock widget for selecting options.
        create_new_job_widget: Creates the widget that allows the user to select a new job.
        wrap_with_border: Wraps a given widget with a border layout.
        _on_item_clicked: Handles the event when an item in the icon grid is clicked.
        clear_dock_widget: Clears the contents of the current dock widget.
        show_job_status_page: Displays the job status page in a dock widget.
        show_new_job_page: Displays the new job selection page in a dock widget.
        clear_layout: Clears the contents of a given layout.
    """
    def __init__(self, parent, job_manager):
        """
        Initialize the MainWidget with a parent widget and job manager.

        Args:
            parent: The parent widget that contains the MainWidget.
            job_manager: The manager responsible for handling job execution.
        """
        super().__init__(parent)
        self.parent = parent
        self.job_manager = job_manager

        # Main layout
        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop | Qt.AlignHCenter)

        # Create the welcome page
        self.clear_layout(self.layout)
        self.new_job_dock = None
        self.dock_widgets = {}

        # Set the main layout
        self.setLayout(self.layout)
 
    def create_dock_widget(self):
        """
        Creates a new dock widget for selecting options. If an existing dock widget 
        is present, it is cleared before creating the new one.
        """
        if self.new_job_dock:
            self.clear_dock_widget()
        self.new_job_dock = QDockWidget("Select options", self)
        self.new_job_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.RightDockWidgetArea)
        self.new_job_dock.setWidget(self.create_new_job_widget())

        if self.parent:
            self.parent.addDockWidget(Qt.RightDockWidgetArea, self.new_job_dock)
    
    def create_new_job_widget(self):
        """
        Creates a widget that allows the user to select a new job.

        Returns:
            new_job_widget: A QWidget that contains the options for selecting a new job.
        """
        new_job_widget = QWidget()
        # Setting specific style for this widget
        new_job_widget.setStyleSheet("""
            QWidget {
                background-color: #CCCCCC;
            }
        """)
        layout = QVBoxLayout(new_job_widget)
        widget = QLabel("Select:")
        font = widget.font()
        font.setPointSize(14)
        widget.setFont(font)
        widget.setMinimumWidth(800)
        widget.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        layout.addWidget(widget)
        self.icon_grid = Icons(self)
        self.icon_grid.addItems(CATEGORIES)
        self.icon_grid.itemClicked.connect(self._on_item_clicked)
        new_job_widget.setMinimumSize(600, 400)

        layout.addWidget(self.icon_grid)
        return new_job_widget
    
    def wrap_with_border(self, widget):
        """
        Wraps a given widget with a border layout.

        Args:
            widget: The QWidget to be wrapped.

        Returns:
            container_widget: A QWidget that wraps the input widget in a border layout.
        """
        container_widget = QWidget()
        layout = QVBoxLayout(container_widget)
        layout.addWidget(widget)
        return container_widget
    
    def _on_item_clicked(self, item):
        """
        Handles the event when an item in the icon grid is clicked. It creates 
        a dock widget with the selected category's options.

        Args:
            item: The item clicked in the icon grid.
        """
        if self.new_job_dock:
            self.new_job_dock.setVisible(False)
        self.create_dock_widget()
        name = item.text()
        widget = CATEGORIES[name].widget(self.job_manager)
        widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.new_job_dock.widget().layout().addWidget(self.wrap_with_border(widget))

    def clear_dock_widget(self):
        """
        Clears the contents of the current dock widget, removing all its child widgets.
        """
        dock_widget = self.new_job_dock.widget()
        layout = dock_widget.layout()
        for i in reversed(range(layout.count())):
            item = layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.deleteLater()
    
    def show_job_status_page(self):
        """
        Displays the job status page in a dock widget.
        """
        job_status_page = JobStatusPage(self.job_manager)
        self.parent.show_dock_widget("Job Status", lambda: job_status_page)
    
    def show_new_job_page(self):
        """
        Displays the new job selection page in a dock widget. If the dock widget 
        already exists, it is made visible again.
        """
        if self.parent:
            self.parent.toolbar.setVisible(True)
        if self.new_job_dock:
            self.new_job_dock.setVisible(False)
        self.create_dock_widget()

    def clear_layout(self, layout, start_index=0):
        """
        Clears the contents of a given layout, starting from a specified index.

        Args:
            layout: The QLayout to clear.
            start_index: The index from which to start clearing the layout.
        """
        while layout.count() > start_index:
            item = layout.takeAt(start_index)
            widget = item.widget()
            if widget:
                widget.deleteLater()
