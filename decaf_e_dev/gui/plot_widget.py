import sys
import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget
)
from pyqtgraph.colormap import ColorMap
from matplotlib.colors import Normalize

class PlotWidget(pg.GraphicsLayoutWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plots = []
        self.setBackground('w')

    def add_plot(self, x_data, y_data, title, x_label, y_label, color=None, label=None, resids=None, scatter=False, colorbar=False):
        plot_item = self.addPlot(title=title)
        plot_item.addLegend()
        self.setBackground('w')
        if scatter and colorbar:
            scatter = self.add_scatter(plot_item, x_data, y_data, resids, colorbar=True)
            colorbar = self.add_colorbar(scatter, resids)
        elif scatter:
            scatter = self.add_scatter(plot_item, x_data, y_data, resids)
        else:
            self.add_line(plot_item, x_data, y_data, color, label)
        
        plot_item.setLabel('left', y_label)
        plot_item.setLabel('bottom', x_label)
        plot_item.showGrid(x=True, y=True)
        self.nextRow()
        self.plots.append(plot_item)
        
        return plot_item
    

    def add_scatter(self, plot_item, x_data, y_data, resids, colorbar=False):
        if colorbar:
            colors = np.array([[68, 1, 84, 255], [58, 82, 139, 255], [32, 144, 140, 255], [94, 201, 97, 255], [253, 231, 37, 255]])
            norm = Normalize(vmin=resids.min(), vmax=resids.max())
            colormap = ColorMap(pos=np.linspace(0, 1, len(colors)), color=colors)
            brushes = [pg.mkBrush(*colormap.map(norm(resid))) for resid in resids]
            
        scatter = pg.ScatterPlotItem(x=x_data, y=y_data, symbol='o', symbolSize=5, brush=brushes, pen=None)
        plot_item.addItem(scatter)
        
        return scatter
    
    def add_line(self, plot_item, x_data, y_data, color, label, lstyle=None):
        plot_item.plot(x_data, y_data, pen=pg.mkPen(color=color, width=2), name=label, linestyle=lstyle)

    def add_colorbar(self, scatter, resids):
        colors = np.array([[68, 1, 84, 255], [58, 82, 139, 255], [32, 144, 140, 255], [94, 201, 97, 255], [253, 231, 37, 255]])
        positions = np.linspace(0, 1, len(colors))
        colormap = ColorMap(pos=positions, color=colors)
        
        colorbar = pg.ColorBarItem(values=(resids.min(), resids.max()), colorMap=colormap, label='Residue #')
        self.addItem(colorbar, row=len(self.plots), col=1)
        
        return colorbar
