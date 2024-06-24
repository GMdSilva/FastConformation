import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from dock_widget import MainWidget
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QStackedWidget, QListWidget, QListWidgetItem, QHBoxLayout, QSizePolicy
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon

class MainFrame(QMainWindow):
    def __init__(self):
        super().__init__()

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        self.layout = QHBoxLayout(self.central_widget)

        # Add an expanding placeholder widget to push the main_widget to the right
        self.placeholder = QWidget()
        self.placeholder.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.layout.addWidget(self.placeholder)

        self.main_widget = MainWidget()
        self.layout.addWidget(self.main_widget)

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

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_frame = MainFrame()
    main_frame.show()
    sys.exit(app.exec_())
