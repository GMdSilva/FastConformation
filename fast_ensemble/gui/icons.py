from pathlib import Path
from qtpy.QtCore import QSize
from qtpy.QtGui import QIcon
from qtpy.QtWidgets import QListWidget, QListWidgetItem
from pathlib import Path
from PyQt5.QtWidgets import QListWidget, QListWidgetItem
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QSize

ICON_ROOT = Path(__file__).parent / "icons"

STYLES = r"""
    QListWidget {
        min-width: 294px;
        background: none;
        font-size: 16pt;
        color: palette(text);
    }
    QListWidget::item {
        width: 95px;
        height: 100px;
        border-radius: 4px;
        margin: 5px;
        padding: 5px;
        background: #D2E3A4;
        color: palette(text);
    }
    QListWidget::item:hover {
        background: palette(dark);
    }
    QListWidget::item:selected {
        background: #ABD149;
        color: palette(highlightedText);
    }
"""

def _get_icon(name):
    """
    Retrieves the file path for an icon based on the given name.

    Args:
        name: The name of the icon file (without extension).

    Returns:
        The file path of the icon as a string if it exists, otherwise an empty string.
    """
    path = ICON_ROOT / f'{name.lower().replace(" ", "_")}.png'
    if not path.exists():
        return ""
    return str(path)

class Icons(QListWidget):
    """
    Icons is a custom QListWidget that displays a grid of icons with labels. 
    It supports adding individual items with optional tooltips and styling.

    Methods:
        addItem: Adds an item with a label and an optional tooltip.
        addItems: Adds multiple items to the list using a dictionary of labels.
    """
    def __init__(self, parent=None):
        """
        Initializes the Icons widget.

        Args:
            parent: The parent widget for this QListWidget.
        """
        super().__init__(parent=parent)
        self.setMovement(self.Static)  # The items cannot be moved by the user.
        self.setViewMode(self.IconMode)  # Make items appear as icons.
        self.setResizeMode(self.Adjust)  # Re-layout when the view is resized.
        self.setUniformItemSizes(True)  # Better performance.
        self.setIconSize(QSize(80, 80))
        self.setWordWrap(True)
        self.setStyleSheet(STYLES)

    def addItem(self, label: str, tool_tip: str = None):
        """
        Adds an item to the QListWidget with a label and an optional tooltip.

        Args:
            label: The label for the item.
            tool_tip: Optional tooltip text for the item.
        """
        if isinstance(label, QListWidgetItem):
            super().addItem(label)
            return

        item = QListWidgetItem(QIcon(_get_icon(label)), label)
        if tool_tip is not None:
            item.setToolTip(tool_tip)
        super().addItem(item)

    def addItems(self, labels: dict) -> None:
        """
        Adds multiple items to the QListWidget from a dictionary of labels.

        Args:
            labels: A dictionary where the keys are labels and the values can contain optional tooltips.
        """
        for label in labels:
            if hasattr(labels[label], "tool_tip"):
                self.addItem(label, labels[label].tool_tip)
            else:
                self.addItem(label)
