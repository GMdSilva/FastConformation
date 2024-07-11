import sys
import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QMessageBox, QScrollArea, QListWidgetItem, QHBoxLayout, QSizePolicy, QPushButton, QToolBar, QDockWidget
)
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Interactive Plot with pyqtgraph')
        self.setGeometry(100, 100, 800, 600)
        self.initUI()

    def initUI(self):
        centralWidget = QWidget()
        layout = QVBoxLayout()
        centralWidget.setLayout(layout)
        self.setCentralWidget(centralWidget)
        imv = pg.ImageView()
        # Create a PlotWidget
        self.plotWidget = pg.PlotWidget()
        layout.addWidget(self.plotWidget)
        
        # Sample data to mimic the provided plot
        x = np.random.uniform(0, 13, 100)
        y = np.random.uniform(40, 100, 100)
        c = np.linspace(200, 450, 100)
        
        # Create scatter plot
        scatter = pg.ScatterPlotItem(x=x, y=y, size=10, pen=None,
                                     brush=pg.mkBrush(255, 255, 255, 120))
        self.plotWidget.addItem(scatter)
        
        # Configure plot axes
        self.plotWidget.setTitle('abl_wt 64 128')
        self.plotWidget.setLabel('left', 'Average pLDDT')
        self.plotWidget.setLabel('bottom', 'C-Alpha RMSF (A)')
        
         
        # Create a colormap
        colors = np.array([[68, 1, 84, 255], [58, 82, 139, 255], [32, 144, 140, 255], [94, 201, 97, 255], [253, 231, 37, 255]])
        positions = np.linspace(0, 1, len(colors))
        colormap = pg.ColorMap(pos=positions, color=colors)
        brushes = [pg.mkBrush(*colormap.map(c_val, mode='float')) for c_val in c]
        scatter.setBrush(brushes)

        # Add color bar
        cbar = pg.ColorBarItem(values=(200, 450), colorMap=colormap)
        # Add interactive annotations
        scatter.sigClicked.connect(self.onScatterClicked)

    def onScatterClicked(self, plot, points):
        for point in points:
            x, y = point.pos()
            index = np.where((plot.data['x'] == x) & (plot.data['y'] == y))[0][0]
            residue = plot.data['data'][index]
            QMessageBox.information(self, 'Point Information', f'Index: {index}\nRMSF: {x:.2f}\npLDDT: {y:.2f}\nResidue #: {residue:.2f}')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWin = MainWindow()
    mainWin.show()
    sys.exit(app.exec_())
