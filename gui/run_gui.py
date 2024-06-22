import sys
from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from dock_widget import MainWidget
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QStackedWidget, QListWidget, QListWidgetItem, QHBoxLayout
)
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtGui import QIcon


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
