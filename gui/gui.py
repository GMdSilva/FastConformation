import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QStackedWidget, QListWidget, QListWidgetItem, QHBoxLayout
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon

@dataclass
class Category:
    widget: Callable
    tool_tip: str = ""

# Placeholder widget classes
class DirectorySelector(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        label = QLabel("Directory Selector")
        layout.addWidget(label)
        self.setLayout(layout)

class PreProcessing(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        label = QLabel("Pre-Processing")
        layout.addWidget(label)
        self.setLayout(layout)

# Define other placeholder widgets as needed
CATEGORIES = {
    "Select files": Category(widget=DirectorySelector, tool_tip="Generating mip and heatmaps"),
    "Pre-Processing": Category(widget=PreProcessing, tool_tip="Preprocessing"),
    # Add other categories here
}

ICON_ROOT = Path(__file__).parent / "icons"
STYLES = r"""
    QListWidget{
        min-width: 100px;
        background: none;
        font-size: 8pt;
        color: #eee;
    }
    QListWidget::item {
        width: 68px;
        height: 100px;
        border-radius: 0;
        margin: 1px;
        padding: 4px;
        background: #414851;
    }
    QListWidget::item::hover {
        background: #5A626C;
    }
"""

def _get_icon(name):
    path = ICON_ROOT / f'{name.lower().replace(" ", "_")}.png'
    if not path.exists():
        return ""
    return str(path)

class Icons(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.setMovement(QListWidget.Static)  # The items cannot be moved by the user.
        self.setViewMode(QListWidget.IconMode)  # make items icons
        self.setResizeMode(QListWidget.Adjust)  # re-layout when view is resized.
        self.setUniformItemSizes(True)  # better performance
        self.setIconSize(QSize(64, 64))
        self.setWordWrap(True)
        self.setStyleSheet(STYLES)

    def addItem(self, label: str, tool_tip: str = None):
        if isinstance(label, QListWidgetItem):
            super().addItem(label)

        item = QListWidgetItem(QIcon(_get_icon(label)), label)
        if tool_tip is not None:
            item.setToolTip(tool_tip)
        super().addItem(item)

    def addItems(self, labels: dict) -> None:
        for label in labels:
            if hasattr(labels[label], "tool_tip"):
                self.addItem(label, labels[label].tool_tip)
            else:
                self.addItem(label)

class MainWidget(QWidget):
    def __init__(self):
        super().__init__()
        layout = QHBoxLayout()
        
        # Sidebar for icons
        sidebar_layout = QVBoxLayout()
        sidebar_layout.setAlignment(Qt.AlignTop)

        widget = QLabel("Select")
        font = widget.font()
        font.setPointSize(10)
        widget.setFont(font)
        widget.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        sidebar_layout.addWidget(widget)

        self.icon_grid = Icons(self)
        self.icon_grid.addItems(CATEGORIES)
        self.icon_grid.itemClicked.connect(self._on_item_clicked)
        sidebar_layout.addWidget(self.icon_grid)

        sidebar_widget = QWidget()
        sidebar_widget.setLayout(sidebar_layout)
        sidebar_widget.setFixedWidth(120)
        
        layout.addWidget(sidebar_widget)

        # Main content area for displaying selected widgets
        self.stacked_widget = QStackedWidget()
        layout.addWidget(self.stacked_widget)

        self.setLayout(layout)

    def _on_item_clicked(self, item):
        name = item.text()
        widget_class = CATEGORIES[name].widget

        for i in range(self.stacked_widget.count()):
            widget = self.stacked_widget.widget(i)
            if type(widget) == widget_class:
                self.stacked_widget.setCurrentWidget(widget)
                return

        new_widget = widget_class()
        self.stacked_widget.addWidget(new_widget)
        self.stacked_widget.setCurrentWidget(new_widget)

class MainFrame(QMainWindow):
    def __init__(self):
        super().__init__()

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout(self.central_widget)
        self.main_widget = MainWidget()
        self.layout.addWidget(self.main_widget)

        self.setWindowTitle('DECAF_E')
        self.showFullScreen()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_frame = MainFrame()
    main_frame.show()
    sys.exit(app.exec_())
