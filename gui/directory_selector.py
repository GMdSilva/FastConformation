import os
from qtpy.QtWidgets import QLabel, QVBoxLayout, QWidget, QApplication, QFileDialog, QListWidget, QPushButton

class DirectorySelector(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        self.directory_label = QPushButton("Select Directory")
        self.directory_label.clicked.connect(self.select_directory)

        self.file_list_label = QLabel("Selected Files:")
        self.file_list_widget = QListWidget()
        self.file_list_widget.doubleClicked.connect(self.add_to_viewer)

        layout.addWidget(self.directory_label)
        layout.addWidget(self.file_list_label)
        layout.addWidget(self.file_list_widget)

        self.setLayout(layout)
        self.setWindowTitle("Directory Selector")

        self.selected_directory = ""

    def select_directory(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ShowDirsOnly | QFileDialog.DontUseNativeDialog
        directory = QFileDialog.getExistingDirectory(
            self, "Select Directory", options=options
        )

        if directory:
            self.selected_directory = directory
            self.directory_label.setText(directory)

            # Read and display files in the selected directory
            self.file_list_widget.clear()
            for file_name in os.listdir(directory):
                if file_name.lower().endswith((".tif", ".tiff", ".nd2")):
                    self.file_list_widget.addItem(file_name)

    def add_to_viewer(self):
        selected_item = self.file_list_widget.currentItem()
        if selected_item:
            file_name = selected_item.text()
            file_path = os.path.join(self.selected_directory, file_name)
            # Here, you need to define your own read function
            # image_data = read(file_path)
            # For this example, I'll just print the file path
            print(f"Selected file: {file_path}")
            # Add your logic to handle the file here