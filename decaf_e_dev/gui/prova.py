import sys
import zmq
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton
from pymol.Qt import get_pymol_widget
import pymol
import threading

class PyMOLWindow(QWidget):
    def __init__(self, parent=None):
        super(PyMOLWindow, self).__init__(parent)
        self.pymol = pymol.PyMOL()
        self.pymol.start()

        layout = QVBoxLayout(self)
        self.pymol_widget = get_pymol_widget(self.pymol)
        layout.addWidget(self.pymol_widget)
        self.setLayout(layout)

        # Initialize ZeroMQ
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PULL)
        self.socket.bind("tcp://*:5555")

        # Start ZeroMQ listener thread
        self.listener_thread = threading.Thread(target=self.listen_to_core)
        self.listener_thread.start()

    def listen_to_core(self):
        while True:
            message = self.socket.recv_string()
            # Process the message from the PySSA core
            self.process_message(message)

    def process_message(self, message):
        # Here you would implement the processing of messages
        print(f"Received message: {message}")

    def load_structure(self, file_path):
        self.pymol.cmd.load(file_path)

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle("PyMOL Integration with PyQt and ZeroMQ")

        # Create the central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create and add the PyMOL window
        self.pymol_window = PyMOLWindow(self)
        layout.addWidget(self.pymol_window)

        # Add a button to send a message to the auxiliary PyMOL instances
        self.send_button = QPushButton("Send Message to Auxiliary Instances")
        self.send_button.clicked.connect(self.send_message)
        layout.addWidget(self.send_button)

        # Initialize ZeroMQ for sending messages
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.PUSH)
        self.socket.connect("tcp://localhost:5556")

    def send_message(self):
        message = "Hello from main PyMOL instance!"
        self.socket.send_string(message)
        print(f"Sent message: {message}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
