from PyQt5.QtWidgets import QWidget, QFormLayout, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QCheckBox, QMessageBox
from PyQt5.QtCore import Qt

def merge_configs(general_options, specific_options):
    config = general_options.copy()
    config.update(specific_options)
    return config

class AnalysisWidgetBase(QWidget):
    def __init__(self):
        super().__init__()


    def validate_inputs(self):
        """Override this method to implement specific validation."""
        return []

    def run_analysis(self, general_options=True):
        errors = self.validate_inputs()
        if errors:
            self.show_error_message(errors)
            return
        
        try:
            g_options = self.general_options_getter()
        except Exception:
            g_options = {}
        
        specific_options = self.get_specific_options()

        config = merge_configs(g_options, specific_options)
        
        try:
            self.run_specific_analysis(config)
        except Exception as e:
            self.show_error_message([str(e)])

    def get_specific_options(self):
        """Override this method to return specific options."""
        return {}

    def run_specific_analysis(self, config):
        """Override this method to run specific analysis."""
        pass

    def show_error_message(self, errors):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Invalid input")
        msg.setInformativeText("\n".join(errors))
        msg.setWindowTitle("Error")
        msg.exec_()
