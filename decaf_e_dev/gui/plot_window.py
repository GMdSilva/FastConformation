from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class PlotWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Plot")
        self.setGeometry(100, 100, 800, 600)
        
        self.main_widget = QWidget(self)
        self.layout = QVBoxLayout(self.main_widget)
        
        self.canvas = FigureCanvas(Figure(figsize=(8, 6)))
        self.layout.addWidget(self.canvas)
        
        self.setCentralWidget(self.main_widget)
        
    def plot(self, plot_func):
        self.canvas.figure.clear()
        ax = self.canvas.figure.subplots()
        plot_func(ax)
        self.canvas.draw()
        self.show()
