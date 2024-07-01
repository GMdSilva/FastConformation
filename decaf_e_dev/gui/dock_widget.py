import warnings
from dataclasses import dataclass
from typing import Callable
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QLabel, QVBoxLayout, QWidget, QApplication, QFileDialog, QListWidget, QPushButton, QHBoxLayout
)
from qtpy.QtGui import QPixmap

# Assuming these imports are correct
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
    def __init__(self):
        super().__init__()

        # Main layout
        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)

        # Create the welcome page
        self.create_welcome_page()

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
        button_layout = QHBoxLayout()
        new_job_button = QPushButton("Submit New Job")
        job_status_button = QPushButton("Job Status")

        new_job_button.clicked.connect(self.show_new_job_page)
        job_status_button.clicked.connect(self.show_job_status_page)

        button_layout.addWidget(new_job_button)
        button_layout.addWidget(job_status_button)
        self.layout.addLayout(button_layout)

    def show_new_job_page(self):
        # Clear the layout
        self.clear_layout(self.layout)

        widget = QLabel("Select")
        font = widget.font()
        font.setPointSize(10)
        widget.setFont(font)
        widget.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.layout.addWidget(widget)

        icon_grid = Icons(self)
        icon_grid.addItems(CATEGORIES)
        icon_grid.itemClicked.connect(self._on_item_clicked)
        self.layout.addWidget(icon_grid)

    def show_job_status_page(self):
        # Placeholder for Job Status page
        self.clear_layout(self.layout)
        label = QLabel("Job Status Page")
        label.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(label)

    def _on_item_clicked(self, item):
        name = item.text()
        widget = CATEGORIES[name].widget()
        font = widget.font()
        font.setPointSize(10)

        # Clear any existing widgets before adding a new one
        self.clear_layout(self.layout, start_index=2)

        # Add the widget to the layout
        self.layout.addWidget(widget)

    def clear_layout(self, layout, start_index=0):
        while layout.count() > start_index:
            item = layout.takeAt(start_index)
            widget = item.widget()
            if widget:
                widget.deleteLater()

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    main_widget = MainWidget()
    main_widget.show()
    sys.exit(app.exec_())
