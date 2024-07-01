import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from decaf_e_dev.gui.dock_widget import MainWidget
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QStackedWidget, QListWidget, QListWidgetItem, QHBoxLayout, QSizePolicy, QPushButton, QToolBar, QAction
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon

class MainFrame(QMainWindow):
    def __init__(self):
        super().__init__()

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QVBoxLayout(self.central_widget)

        # Add the main widget
        self.main_widget = MainWidget(self)
        self.layout.addWidget(self.main_widget)

        # Add the toolbar
        self.toolbar = QToolBar()
        self.toolbar.setMovable(False)
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

        self.home_button.clicked.connect(self.show_home_page)
        self.submit_new_job_button.clicked.connect(self.show_new_job_page)
        self.job_status_button.clicked.connect(self.show_job_status_page)

        self.toolbar.addWidget(self.home_button)
        self.toolbar.addWidget(self.submit_new_job_button)
        self.toolbar.addWidget(self.job_status_button)
        self.toolbar.setVisible(False)  # Initially hidden

        self.setWindowTitle('DECAF_E')
        self.showFullScreen()

        self.apply_styles()

    def apply_styles(self):
        self.setStyleSheet("""
            QWidget {
                background-color: #f5f5f5;
                color: #333333;
            }
            QLabel {
                color: #333333;
            }
            QPushButton {
                background-color: #007BFF;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 5px;
            }
            QPushButton:hover {
                background-color: #0056b3;
            }
            QPushButton:pressed {
                background-color: #004494;
            }
            QLineEdit {
                background-color: white;
                border: 1px solid #CCCCCC;
                padding: 4px;
                border-radius: 4px;
                color: #333333;
            }
            QComboBox {
                background-color: white;
                border: 1px solid #CCCCCC;
                padding: 4px;
                border-radius: 4px;
                color: #333333;
            }
            QListWidget {
                background-color: white;
                border: 1px solid #CCCCCC;
                padding: 4px;
                color: #333333;
            }
        """)

    def show_home_page(self):
        self.main_widget.create_welcome_page()
        self.toolbar.setVisible(False)

    def show_new_job_page(self):
        self.main_widget.show_new_job_page()
        self.toolbar.setVisible(True)

    def show_job_status_page(self):
        self.main_widget.show_job_status_page()
        self.toolbar.setVisible(True)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_frame = MainFrame()
    main_frame.show()
    sys.exit(app.exec_())
