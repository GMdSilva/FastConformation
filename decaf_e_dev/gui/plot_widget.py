import sys
import os
import glob
import pandas as pd
import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QWidget, QFormLayout, QLineEdit, QPushButton, QVBoxLayout, QHBoxLayout, QFileDialog, QApplication
)
from PyQt5.QtCore import Qt
from tqdm import tqdm

class PlotWidget(pg.GraphicsLayoutWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plots = []

    def add_plot(self, x_data, y_data, title, x_label, y_label):
        plot_item = self.addPlot(title=title)
        plot_item.plot(x_data, y_data, pen=pg.mkPen(color='b', width=2))
        plot_item.setLabel('left', y_label)
        plot_item.setLabel('bottom', x_label)
        plot_item.showGrid(x=True, y=True)
        self.nextRow()
        self.plots.append(plot_item)
