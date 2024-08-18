import sys
import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget
)
from pyqtgraph.colormap import ColorMap
from matplotlib.colors import Normalize

class PlotWidget(pg.GraphicsLayoutWidget):
    """
    PlotWidget is a custom widget for plotting data using pyqtgraph. It supports
    both line plots and scatter plots with optional color mapping and colorbars.

    Methods:
        add_plot: Adds a new plot to the widget.
        add_borders: Adds borders to the plot.
        add_scatter: Adds a scatter plot to the plot item.
        add_line: Adds a line plot to the plot item.
        add_colorbar: Adds a colorbar to the widget.
    """
    def __init__(self, parent=None):
        """
        Initializes the PlotWidget.

        Args:
            parent: The parent widget, if any.
        """
        super().__init__(parent)
        self.plots = []
        self.setBackground('w')

    def add_plot(self, x_data, y_data, title, x_label, y_label, color=None, label=None, resids=None, scatter=False, colorbar=False):
        """
        Adds a new plot to the widget.

        Args:
            x_data: The data for the x-axis.
            y_data: The data for the y-axis.
            title: The title of the plot.
            x_label: The label for the x-axis.
            y_label: The label for the y-axis.
            color: The color of the line or scatter points.
            label: The label for the legend.
            resids: Residual values for coloring scatter points.
            scatter: Whether to create a scatter plot.
            colorbar: Whether to add a colorbar for the scatter plot.

        Returns:
            The created plot item.
        """
        plot_item = self.addPlot(title=title)
        plot_item.addLegend()
        self.setBackground('w')
        if scatter and colorbar:
            scatter = self.add_scatter(plot_item, x_data, y_data, resids, colorbar=True)
            colorbar = self.add_colorbar(resids)
        elif scatter:
            scatter = self.add_scatter(plot_item, x_data, y_data, resids)
        else:
            self.add_line(plot_item, x_data, y_data, color, label)
        
        plot_item.setLabel('left', y_label)
        plot_item.setLabel('bottom', x_label)
        plot_item.showGrid(x=True, y=True)
        self.nextRow()
        self.plots.append(plot_item)
        self.add_borders(plot_item)
        return plot_item
    
    def add_borders(self, plot):
        """
        Adds borders to the given plot.

        Args:
            plot: The plot item to which borders will be added.
        """
        plot.getViewBox().setBorder(pg.mkPen(color='lightgrey', width=1))

    def add_scatter(self, plot_item, x_data, y_data, resids=None, color='b', colorbar=False, label=None):
        """
        Adds a scatter plot to the given plot item.

        Args:
            plot_item: The plot item to which the scatter plot will be added.
            x_data: The data for the x-axis.
            y_data: The data for the y-axis.
            resids: Residual values for coloring the scatter points.
            color: The color of the scatter points.
            colorbar: Whether to add a colorbar.
            label: The label for the legend.

        Returns:
            The created scatter plot item.
        """
        if colorbar and (resids is not None):
            colors = np.array([[68, 1, 84, 255], [58, 82, 139, 255], [32, 144, 140, 255], [94, 201, 97, 255], [253, 231, 37, 255]])
            norm = Normalize(vmin=resids.min(), vmax=resids.max())
            colormap = ColorMap(pos=np.linspace(0, 1, len(colors)), color=colors)
            brushes = [pg.mkBrush(*colormap.map(norm(resid))) for resid in resids]
        else:
            pg_color = pg.mkColor(color)
            brushes = pg.mkBrush(pg_color)
        scatter = pg.ScatterPlotItem(x=x_data, y=y_data, symbol='o', symbolSize=5, brush=brushes, name=label, pen=None)
        plot_item.addItem(scatter)
        
        return scatter
    
    def add_line(self, plot_item, x_data, y_data, color, label, lstyle=None):
        """
        Adds a line plot to the given plot item.

        Args:
            plot_item: The plot item to which the line plot will be added.
            x_data: The data for the x-axis.
            y_data: The data for the y-axis.
            color: The color of the line.
            label: The label for the legend.
            lstyle: The line style (e.g., solid, dashed).
        """
        plot_item.plot(x_data, y_data, pen=pg.mkPen(color=color, width=2), name=label, linestyle=lstyle)

    def add_colorbar(self, resids):
        """
        Adds a colorbar to the widget based on residual values.

        Args:
            resids: Residual values to map the colorbar.

        Returns:
            The created colorbar item.
        """
        colors = np.array([[68, 1, 84, 255], [58, 82, 139, 255], [32, 144, 140, 255], [94, 201, 97, 255], [253, 231, 37, 255]])
        positions = np.linspace(0, 1, len(colors))
        colormap = ColorMap(pos=positions, color=colors)
        
        colorbar = pg.ColorBarItem(values=(resids.min(), resids.max()), colorMap=colormap, label='Residue #')
        self.addItem(colorbar, row=len(self.plots), col=1)
        
        return colorbar