import warnings
from dataclasses import dataclass
from typing import Callable
from decaf_e_dev.gui.directory_selector import DirectorySelector
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QLabel, QVBoxLayout, QWidget, QApplication, QFileDialog, QListWidget, QPushButton
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
        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignTop)

        widget = QLabel("Select")
        font = widget.font()
        font.setPointSize(10)
        widget.setFont(font)
        widget.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        layout.addWidget(widget)

        icon_grid = Icons(self)
        icon_grid.addItems(CATEGORIES)
        icon_grid.itemClicked.connect(self._on_item_clicked)
        layout.addWidget(icon_grid)

        self.setLayout(layout)
        self.setWindowTitle("DECAF")
    
    def _on_item_clicked(self, item):
        name = item.text()
        widget = CATEGORIES[name].widget()
        font = widget.font()
        font.setPointSize(10)
        
        # Clear any existing widgets before adding new one
        while self.layout().count() > 2:
            item = self.layout().takeAt(2)
            item.widget().deleteLater()

        # Add the widget to the layout
        layout = self.layout()
        layout.addWidget(widget)