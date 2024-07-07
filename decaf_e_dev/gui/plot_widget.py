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
from pyqtgraph.colormap import ColorMap
class PlotWidget(pg.GraphicsLayoutWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plots = []
        self.setBackground('w')

    def add_plot(self, x_data, y_data, title, x_label, y_label, color=None, label=None, cmap=None, resids=None, scatter=False):
        
        plot_item = self.addPlot(title=title)
        plot_item.addLegend()
        self.setBackground('w')
        plot_item.setFixedHeight(600)
        plot_item.setFixedWidth(800)
        if not (cmap and resids) is None:
            colorbar=self.add_colorbar(resids, cmap) 
        if scatter:
            scatter=self.add_scatter(plot_item, x_data, y_data, color, label)
            if not (cmap and resids) is None:
                colorbar.setImageItem(scatter)
        else:
            self.add_line(plot_item, x_data, y_data, color, label)
        
        plot_item.setLabel('left', y_label)
        plot_item.setLabel('bottom', x_label)
        plot_item.showGrid(x=True, y=True)
        self.nextRow()
        self.plots.append(plot_item)
        
        return plot_item
    
    def add_line(self, plot_item, x_data, y_data, color, label):
        plot_item.plot(x_data, y_data, pen=pg.mkPen(color=color, width=2), name=label)

    def add_scatter(self, plot_item, x_data, y_data, color_data, label=None):
        scatter = pg.ScatterPlotItem(x=x_data, y=y_data, pen=None, brush=color_data, name=label)
        plot_item.addItem(scatter)
        return scatter

    def add_colorbar(self, resids, cmap):
        # Create colorbar
        colormap = ColorMap(pos=np.linspace(0, 1, cmap.N), color=cmap.colors)
        colorbar = pg.ColorBarItem(values=(resids.min(), resids.max()), colorMap=colormap, interactive=True)
        self.addItem(colorbar, row=1, col=1)
        return colorbar
