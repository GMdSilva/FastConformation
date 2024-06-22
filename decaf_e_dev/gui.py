import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QFormLayout, QLineEdit, QPushButton, QFileDialog, QLabel
from PyQt5.QtGui import QPixmap
import json

class MainWindow(QMainWindow):
    def __init__(self, config):
        super().__init__()

        self.config = config
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Config GUI')
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout(central_widget)

        self.form_layout = QFormLayout()
        self.fields = {}

        for key, value in self.config.items():
            if key.startswith('_comment'):
                continue
            
            if isinstance(value, list):
                value = ','.join(map(str, value))
            elif isinstance(value, bool):
                value = str(value)

            label = QLabel(key)
            line_edit = QLineEdit(str(value))
            self.fields[key] = line_edit
            self.form_layout.addRow(label, line_edit)

        layout.addLayout(self.form_layout)

        save_button = QPushButton('Save Config')
        save_button.clicked.connect(self.save_config)
        layout.addWidget(save_button)

        self.filepath_button = QPushButton('Select File')
        self.filepath_button.clicked.connect(self.open_file_dialog)
        layout.addWidget(self.filepath_button)

        self.image_label = QLabel(self)
        layout.addWidget(self.image_label)

    def save_config(self):
        new_config = {}
        for key, line_edit in self.fields.items():
            value = line_edit.text()
            if ',' in value:
                value = list(map(str.strip, value.split(',')))
            elif value.lower() in ['true', 'false']:
                value = value.lower() == 'true'
            new_config[key] = value

        with open('new_config.json', 'w') as file:
            json.dump(new_config, file, indent=4)

    def open_file_dialog(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", "All Files (*);;PDB Files (*.pdb)", options=options)
        if file_path:
            if file_path.endswith('.pdb'):
                self.display_pdb(file_path)

    def display_pdb(self, file_path):
        # Placeholder for displaying PDB content
        # You can use tools like PyMOL, nglview, or other visualization libraries
        pixmap = QPixmap(file_path)  # Placeholder for actual PDB rendering
        self.image_label.setPixmap(pixmap)
        self.image_label.setScaledContents(True)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    with open('example_configs/config.json') as f:
        config = json.load(f)
    main_win = MainWindow(config)
    main_win.show()
    sys.exit(app.exec_())
