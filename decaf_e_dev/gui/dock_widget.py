from dataclasses import dataclass
from typing import Callable
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QLabel, QVBoxLayout, QWidget, QDockWidget, QPushButton, QHBoxLayout, QSizePolicy
)
from qtpy.QtGui import QPixmap
from decaf_e_dev.gui.directory_selector import DirectorySelector
from decaf_e_dev.gui.icons import Icons
from decaf_e_dev.gui.build_msa import MSAOptionsWidget
from decaf_e_dev.gui.make_predictions import MakePredictionsWidget
from decaf_e_dev.gui.analysis_config import AnalysisConfigWidget

@dataclass
class Category:
    widget: Callable
    tool_tip: str = ""

CATEGORIES = {
    "Build MSA": Category(
        widget=lambda: MSAOptionsWidget(),
        tool_tip="Select parameters to build MSA",
    ),
    "Make Predictions": Category(
        widget=lambda: MakePredictionsWidget(),
        tool_tip="Select parameters to make predictions",
    ),
    "Analysis": Category(
        widget=lambda: AnalysisConfigWidget(),
        tool_tip="Select parameters to analyze results",
    ),
}

class MainWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent

        # Main layout
        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignCenter)

        # Create the welcome page
        self.create_welcome_page()

        self.dock_widgets={}
        
        # Set the main layout
        self.setLayout(self.layout)
        self.setWindowTitle("decaf_dev")


    def create_welcome_page(self):
        # Clear the layout
        self.clear_layout(self.layout)

        # Add the title
        title = QLabel("decaf_dev")
        font = title.font()
        font.setPointSize(20)
        title.setFont(font)
        title.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(title)

        # Add the image
        image = QLabel(self)
        pixmap = QPixmap('/Users/fmgaleazzi/decaf_e_dev/methods-2.png')  # Change this to the path of your image
        image.setPixmap(pixmap)
        image.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(image)

        # Add the buttons
        self.button_layout = QHBoxLayout()
        self.new_job_button = QPushButton("Submit New Job")
        self.job_status_button = QPushButton("Job Status")

        self.new_job_button.clicked.connect(self.show_new_job_page)
        self.job_status_button.clicked.connect(self.show_job_status_page)
        self.new_job_dock=None
        self.button_layout.addWidget(self.new_job_button)
        self.button_layout.addWidget(self.job_status_button)
        self.layout.addLayout(self.button_layout)

    def create_dock_widget(self):
        if self.new_job_dock:
            self.clear_dock_widget()
        self.clear_layout(self.button_layout)
        self.new_job_dock = QDockWidget("Select options", self)
        self.new_job_dock.setAllowedAreas(Qt.TopDockWidgetArea | Qt.RightDockWidgetArea)
        self.new_job_dock.setWidget(self.create_new_job_widget())
        self.new_job_dock.setVisible(True)
        if self.parent:
            self.parent.addDockWidget(Qt.RightDockWidgetArea, self.new_job_dock)
    
    def create_new_job_widget(self):
        new_job_widget = QWidget()
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
    
    def _on_item_clicked(self, item):
        name = item.text()
        widget = CATEGORIES[name].widget()
        widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.new_job_dock.widget().layout().addWidget(widget)

    def clear_dock_widget(self):
        dock_widget = self.new_job_dock.widget()
        layout = dock_widget.layout()
        for i in reversed(range(layout.count())):
            item = layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.deleteLater()
    
    def show_new_job_page(self):
        if self.parent:
            self.parent.toolbar.setVisible(True)
        self.create_dock_widget()

    def show_job_status_page(self):
        self.central_widget.show_job_status_page()
        self.toolbar.setVisible(True)
        self.hide_all_dock_widgets()
    
    def _on_item_clicked(self, item):
        name = item.text()
        widget_class = CATEGORIES[name].widget
        self.parent.show_dock_widget(name, widget_class)
        
    def clear_layout(self, layout, start_index=0):
        while layout.count() > start_index:
            item = layout.takeAt(start_index)
            widget = item.widget()
            if widget:
                widget.deleteLater()
