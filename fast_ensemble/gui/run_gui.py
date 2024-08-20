import sys
from dataclasses import dataclass
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from PyQt5.QtWidgets import (
    QApplication, QDesktopWidget, QMainWindow, QVBoxLayout, QWidget, QPlainTextEdit, QScrollArea, QHBoxLayout, QSizePolicy, QPushButton, QToolBar, QDockWidget, QLabel
)
import sys
sys.path.append('/Users/fmgaleazzi/fast_ensemble')
from fast_ensemble.gui.job_manager import JobStatusPage, JobManager
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout
from PyQt5.QtGui import QPixmap, QPainter, QIcon
from PyQt5.QtCore import Qt
from fast_ensemble.gui.dock_widget import MainWidget
from fast_ensemble.gui.icons import Icons
from fast_ensemble.gui.build_msa import MSAOptionsWidget
from fast_ensemble.gui.make_predictions import MakePredictionsWidget
from fast_ensemble.gui.analysis_config import AnalysisConfigWidget
import signal
import warnings

class MainFrame(QMainWindow):
    """
    MainFrame is the main application window for the FastEnsemble application. 
    It contains the central widget, toolbar, and dock widgets for managing different tasks.

    Methods:
        create_terminal_dock: Creates a dock widget for displaying the analysis log (terminal output).
        show_terminal: Toggles the visibility of the terminal dock widget.
        apply_styles: Applies custom styles to the main window and widgets.
        show_home_page: Displays the home page, hiding other dock widgets.
        show_new_job_page: Displays the new job page, allowing users to submit new jobs.
        show_job_status_page: Displays the job status page, showing the status of submitted jobs.
        show_plot: Displays a plot in a dock widget.
        show_dock_widget: Displays a specified dock widget.
        hide_all_dock_widgets: Hides all currently visible dock widgets.
        set_initial_window_size: Sets the initial size of the main window.
    """
    def __init__(self):
        """
        Initializes the MainFrame with its central widget, toolbar, and job manager.
        """
        super().__init__()
        self.job_manager = JobManager()
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout(self.central_widget)
        self.background_widget = BackgroundWidget(self)
        self.layout.addWidget(self.background_widget)

        # Add the main widget
        self.main_widget = MainWidget(self, self.job_manager)
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
        self.job_status_button = QPushButton("MSA and Predictions Job Status")
        self.terminal_button = QPushButton("Analysis Log")

        self.home_button.clicked.connect(self.show_home_page)
        self.submit_new_job_button.clicked.connect(self.show_new_job_page)
        self.job_status_button.clicked.connect(self.show_job_status_page)
        self.terminal_button.clicked.connect(self.show_terminal)

        self.toolbar.addWidget(self.home_button)
        self.toolbar.addWidget(self.submit_new_job_button)
        self.toolbar.addWidget(self.job_status_button)
        self.toolbar.addWidget(self.terminal_button)
        self.toolbar.setVisible(True)
        
        self.setWindowTitle('FastEnsemble')
        self.set_initial_window_size()

        self.apply_styles()

        # Redirect stdout and stderr
        self.terminal_dock = None
        self.create_terminal_dock()
        sys.stdout = QPlainTextEditLogger(self.terminal_output)
        sys.stderr = QPlainTextEditLogger(self.terminal_output)

    def create_terminal_dock(self):
        """
        Creates a dock widget for displaying the analysis log (terminal output).
        """
        self.terminal_output = QPlainTextEdit(self)
        self.terminal_output.setReadOnly(True)

        self.terminal_dock = QDockWidget("Show Analysis Log", self)
        self.terminal_dock.setAllowedAreas(Qt.LeftDockWidgetArea)
        self.terminal_dock.setWidget(self.terminal_output)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.terminal_dock)
        self.terminal_dock.setVisible(False)

    def show_terminal(self):
        """
        Toggles the visibility of the terminal dock widget.
        """
        if self.terminal_dock.isVisible():
            self.terminal_dock.setVisible(False)
        else:
            self.terminal_dock.setVisible(True)

    def apply_styles(self):
        """
        Applies custom styles to the main window and its widgets.
        """
        self.setStyleSheet("""
            QWidget {
                background-color: palette(base);
                color: palette(text);
                font-size: 16px;
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
                font-size: 16px;
            }
            QComboBox {
                background-color: palette(base);
                border: 1px solid palette(mid);
                padding: 4px;
                border-radius: 4px;
                color: palette(text);
                font-size: 16px;
            }
            QListWidget {
                background-color: palette(base);
                border: 1px solid palette(mid);
                padding: 4px;
                color: palette(text);
                font-size: 16px;
            }
            QDockWidget {
                background-color: darkgrey;
            }
        """)

    def show_home_page(self):
        """
        Displays the home page, hiding all other dock widgets.
        """
        self.toolbar.setVisible(True)
        self.hide_all_dock_widgets()
        self.main_widget.new_job_dock.setVisible(False)
        self.terminal_dock.setVisible(False)

    def show_new_job_page(self):
        """
        Displays the new job submission page, hiding other dock widgets.
        """
        self.hide_all_dock_widgets()
        if self.main_widget.new_job_dock:
            self.main_widget.new_job_dock.setVisible(False)
        self.toolbar.setVisible(True)
        self.main_widget.create_dock_widget()

    def show_job_status_page(self):
        """
        Displays the job status page, hiding other dock widgets.
        """
        self.hide_all_dock_widgets()
        self.main_widget.new_job_dock.setVisible(False)
        self.toolbar.setVisible(True)
        self.main_widget.show_job_status_page()

    def show_plot(self, title, content_widget):
        """
        Displays a plot in a new dock widget.

        Args:
            title: The title of the dock widget.
            content_widget: The widget containing the plot content.
        """
        dock_widget = QDockWidget(title, self)
        dock_widget.setAllowedAreas(Qt.AllDockWidgetAreas)

        container_widget = QWidget()
        layout = QVBoxLayout(container_widget)
        layout.addWidget(content_widget)
        container_widget.setStyleSheet("QWidget { background-color: lightgrey; }")

        dock_widget.setWidget(container_widget)
        dock_widget.raise_()

        self.addDockWidget(Qt.LeftDockWidgetArea, dock_widget)
        dock_widget.setMinimumWidth(800)
        self.dock_widgets[title] = dock_widget

    def show_dock_widget(self, title, widget_callable):
        """
        Displays a specified dock widget. If it does not exist, creates and displays it.

        Args:
            title: The title of the dock widget.
            widget_callable: A callable that returns the content widget for the dock.
        """
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
        """
        Hides all currently visible dock widgets.
        """
        for dock_widget in self.dock_widgets.values():
            dock_widget.setVisible(False)

    def set_initial_window_size(self):
        """
        Sets the initial size of the main window based on the screen size.
        """
        screen = QDesktopWidget().screenGeometry()
        self.setGeometry(100, 100, screen.width() - 100, screen.height() - 100)

class QPlainTextEditLogger:
    """
    A custom logger that redirects stdout and stderr to a QPlainTextEdit widget.

    Methods:
        write: Appends a message to the text edit widget.
        flush: Dummy method to comply with the logging interface.
    """
    def __init__(self, text_edit):
        """
        Initializes the logger with the QPlainTextEdit widget.

        Args:
            text_edit: The QPlainTextEdit widget where logs will be displayed.
        """
        self.widget = text_edit

    def write(self, message):
        """
        Appends a message to the text edit widget.

        Args:
            message: The message to append.
        """
        self.widget.appendPlainText(message)

    def flush(self):
        """
        Dummy method to comply with the logging interface.
        """
        pass

class BackgroundWidget(QWidget):
    """
    A custom QWidget that displays a background image in the main window.

    Methods:
        paintEvent: Handles the painting of the background image.
    """
    def __init__(self, parent=None):
        """
        Initializes the BackgroundWidget with a parent widget.

        Args:
            parent: The parent widget, if any.
        """
        super(BackgroundWidget, self).__init__(parent)
        self.pixmap = QPixmap("background_logo.png")

    def paintEvent(self, event):
        """
        Handles the painting of the background image, centering it in the widget.

        Args:
            event: The paint event.
        """
        painter = QPainter(self)
        pixmap_rect = self.pixmap.rect()
        center_x = (self.width() - pixmap_rect.width()) // 2
        center_y = (self.height() - pixmap_rect.height()) + 30
        painter.drawPixmap(center_x, center_y, self.pixmap)

def handle_sigint(signal, frame):
    """
    Handles the SIGINT signal, allowing the application to shut down gracefully.

    Args:
        signal: The signal received.
        frame: The current stack frame.
    """
    print("SIGINT received, shutting down...")
    QApplication.quit()

def main():
    """
    The main entry point for the application. It sets up the application, handles signals,
    and starts the main event loop.
    """
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    signal.signal(signal.SIGINT, handle_sigint)
    app = QApplication(sys.argv)
    app_icon = QIcon('methods-2.png')  # Update the path as needed
    app.setWindowIcon(app_icon)
    main_frame = MainFrame()
    main_frame.show()
    try:
        sys.exit(app.exec_())
    except KeyboardInterrupt:
        print("Shutting down...")
        main_frame.job_manager.stop()
        sys.exit(0)

if __name__ == '__main__':
    main()