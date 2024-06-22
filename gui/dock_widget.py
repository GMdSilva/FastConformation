import warnings
from dataclasses import dataclass
from typing import Callable
from directory_selector import DirectorySelector
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QLabel, QVBoxLayout, QWidget, QApplication, QFileDialog, QListWidget, QPushButton
from icons import Icons
@dataclass
class Category:
    widget: Callable
    tool_tip: str = ""

# Assuming preprocessing and segmentation are defined elsewhere
# Placeholder for these functions for this example
def preprocessing():
    pass

def segmentation():
    pass

CATEGORIES = {
    "Select Output Path": Category(
        widget=lambda: DirectorySelector(),  # No viewer argument required
        tool_tip="Select an output directory",
    ),
    "Build MSA": Category(
        widget=preprocessing,
        tool_tip="Select parameters to build MSA",
    ),
    "Make Predictions": Category(
        widget=segmentation, 
        tool_tip="Run ensemble prediction"
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
        self.setWindowTitle("Voltage Imaging Analysis")
    
    def _on_item_clicked(self, item):
        name = item.text()
        widget = CATEGORIES[name].widget()
        # Add the widget to the layout
        layout = self.layout()
        layout.addWidget(widget)


if __name__ == "__main__":
    import sys
    import os

    app = QApplication(sys.argv)
    main_window = MainWidget()
    main_window.show()
    sys.exit(app.exec_())