import warnings
from dataclasses import dataclass
from typing import Callable

from qtpy.QtCore import Qt
from qtpy.QtWidgets import QLabel, QVBoxLayout, QWidget

from .directory_selector import DirectorySelector
from .icons import Icons


@dataclass
class Category:
    widget: Callable
    tool_tip: str = ""


CATEGORIES = {
    "Select Output Path": Category(
        widget=DirectorySelector,
        tool_tip="Select an output directory",
    ),
    "Build MSA": Category(
        widget=preprocessing,
        tool_tip="Select parameters to build MSA",
    ),
    "Make Predictions": Category(
        widget=segmentation, tool_tip="Run ensemble prediction"
    ),
}


class MainWidget(QWidget):
    def __init__(self, viewer: napari.viewer.Viewer):
        super().__init__()
        self._viewer = viewer
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
        widget = CATEGORIES[name].widget

        if name == "Plot Neuron Data" or name == "Patch Clamp Visualization":
            try:
                widget = widget(self._viewer)
                self._viewer.window.add_dock_widget(
                    widget, area="bottom", name=name
                )
            except AttributeError:
                napari.utils.notifications.show_warning(
                    "No segmentation layer detected"
                )
        else:
            if name == "Select files":
                widget = widget(self._viewer)
                area = "left"
            else:
                area = "right"
            self._viewer.window.add_dock_widget(widget, area=area, name=name)
