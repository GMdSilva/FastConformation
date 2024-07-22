import sys
from dataclasses import dataclass
import sys
sys.path.append('/Users/fmgaleazzi/decaf_e_dev')
import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from PyQt5.QtWidgets import (
    QApplication, QDesktopWidget, QMainWindow, QVBoxLayout, QWidget, QPlainTextEdit, QScrollArea, QHBoxLayout, QSizePolicy, QPushButton, QToolBar, QDockWidget, QLabel
)
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout
from PyQt5.QtGui import QPixmap, QPainter, QIcon
from PyQt5.QtCore import Qt
from decaf_e_dev.gui.dock_widget import MainWidget
from decaf_e_dev.gui.icons import Icons
from decaf_e_dev.gui.directory_selector import DirectorySelector
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
class BackgroundWidget(QWidget):
    def __init__(self, parent=None):
        super(BackgroundWidget, self).__init__(parent)
        self.pixmap = QPixmap("methods-2.png")

    def paintEvent(self, event):
        painter = QPainter(self)
        pixmap_rect = self.pixmap.rect()
        center_x = (self.width() - pixmap_rect.width()) // 2
        center_y = (self.height() - pixmap_rect.height()) + 50
        painter.drawPixmap(center_x, center_y, self.pixmap)

class MainFrame(QMainWindow):
    def __init__(self):
        super().__init__()

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)


        self.layout = QVBoxLayout(self.central_widget)
        self.background_widget = BackgroundWidget(self)
        self.layout.addWidget(self.background_widget)
        # Add the main widget
        self.main_widget = MainWidget(self)
        self.layout.addWidget(self.main_widget)

        # Add the toolbar
        self.toolbar = QToolBar()
        self.toolbar.setMovable(False)
        self.dock_widgets = {}
        self.toolbar.setStyleSheet("""
            QToolBar {
                background-color: #333333;
                color: white;
                padding: 10px;
            }
            QPushButton {
                background-color: #555555;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 0 5px;
            }
            QPushButton:hover {
                background-color: #666666;
            }
            QPushButton:pressed {
                background-color: #777777;
            }
        """)

        self.addToolBar(Qt.TopToolBarArea, self.toolbar)

        self.home_button = QPushButton("Home Page")
        self.submit_new_job_button = QPushButton("Submit New Job")
        self.job_status_button = QPushButton("Job Status")
        self.terminal_button = QPushButton("Log")

        self.home_button.clicked.connect(self.show_home_page)
        self.submit_new_job_button.clicked.connect(self.show_new_job_page)
        self.job_status_button.clicked.connect(self.show_job_status_page)
        self.terminal_button.clicked.connect(self.show_terminal)

        self.toolbar.addWidget(self.home_button)
        self.toolbar.addWidget(self.submit_new_job_button)
        self.toolbar.addWidget(self.job_status_button)
        self.toolbar.addWidget(self.terminal_button)
        self.toolbar.setVisible(False)

        self.setWindowTitle('DECAF_E')
        self.set_initial_window_size()

        self.apply_styles()

        # Redirect stdout and stderr
        self.terminal_dock = None
        self.create_terminal_dock()
        sys.stdout = QPlainTextEditLogger(self.terminal_output)
        sys.stderr = QPlainTextEditLogger(self.terminal_output)

    def create_terminal_dock(self):
        self.terminal_output = QPlainTextEdit(self)
        self.terminal_output.setReadOnly(True)

        self.terminal_dock = QDockWidget("Show Log", self)
        self.terminal_dock.setAllowedAreas(Qt.LeftDockWidgetArea)
        self.terminal_dock.setWidget(self.terminal_output)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.terminal_dock)
        self.terminal_dock.setVisible(False)

    def show_terminal(self):
        if self.terminal_dock.isVisible():
            self.terminal_dock.setVisible(False)
        else:
            self.terminal_dock.setVisible(True)

    def apply_styles(self):
        self.setStyleSheet("""
            QWidget {
                background-color: palette(base);
                color: palette(text);
            }
            QLabel {
                color: palette(text);
            }
            QPushButton {
                background-color: #D2E3A4;
                color: palette(highlightedText);
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: palette(dark);
            }
            QPushButton:pressed {
                background-color: #ABD149;
            }
            QLineEdit {
                background-color: palette(base);
                border: 1px solid palette(mid);
                padding: 4px;
                border-radius: 4px;
                color: palette(text);
            }
            QComboBox {
                background-color: palette(base);
                border: 1px solid palette(mid);
                padding: 4px;
                border-radius: 4px;
                color: palette(text);
            }
            QListWidget {
                background-color: palette(base);
                border: 1px solid palette(mid);
                padding: 4px;
                color: palette(text);
            }
            QDockWidget {
                background-color: lightgrey;
            }
        """)

    def show_home_page(self):
        self.toolbar.setVisible(False)
        self.hide_all_dock_widgets()
        self.main_widget.new_job_dock.setVisible(False)
        self.main_widget.create_welcome_page()

    def show_new_job_page(self):
        self.hide_all_dock_widgets()
        self.main_widget.new_job_dock.setVisible(False)
        self.toolbar.setVisible(True)
        self.main_widget.create_dock_widget()

    def show_job_status_page(self):
        self.main_widget.new_job_dock.setVisible(False)
        self.main_widget.show_job_status_page()
        self.toolbar.setVisible(True)

    def show_plot(self, title, content_widget):
        # Create a dock widget
        dock_widget = QDockWidget(title, self)
        dock_widget.setAllowedAreas(Qt.AllDockWidgetAreas)

        dock_widget.setWidget(content_widget)
        dock_widget.raise_()
        # Add the dock widget to the main window, ensuring it takes up as much space as possible
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_widget)
        dock_widget.setMinimumWidth(800)
        self.dock_widgets[title] = dock_widget

    def show_dock_widget(self, title, widget_callable):
        if title in self.dock_widgets:
            dock_widget = self.dock_widgets[title]
            dock_widget.setVisible(True)
            dock_widget.raise_()
        else:
            dock_widget = QDockWidget(title, self)
            dock_widget.setAllowedAreas(Qt.RightDockWidgetArea)

            scroll_area = QScrollArea()
            scroll_area.setWidgetResizable(True)
            content_widget = widget_callable()
            scroll_area.setWidget(content_widget)

            dock_widget.setWidget(scroll_area)
            self.addDockWidget(Qt.RightDockWidgetArea, dock_widget)
            self.dock_widgets[title] = dock_widget

        for dock_key, dock_value in self.dock_widgets.items():
            if dock_key != title:
                dock_value.setVisible(False)

    def hide_all_dock_widgets(self):
        for dock_widget in self.dock_widgets.values():
            dock_widget.setVisible(False)

    def set_initial_window_size(self):
        screen = QDesktopWidget().screenGeometry()
        self.setGeometry(100, 100, screen.width() - 100, screen.height() - 100)

# show "terminal" on screen
class QPlainTextEditLogger:
    def __init__(self, text_edit):
        self.widget = text_edit

    def write(self, message):
        self.widget.appendPlainText(message)

    def flush(self):
        pass
    
def main():
    app = QApplication(sys.argv)
    # Set the application icon
    app_icon = QIcon('methods-2.png')  # Update the path as needed
    app.setWindowIcon(app_icon)

    main_frame = MainFrame()
    main_frame.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()