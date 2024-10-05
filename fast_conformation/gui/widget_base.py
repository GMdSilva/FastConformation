from PyQt5.QtWidgets import QWidget, QFormLayout, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QMessageBox
from PyQt5.QtCore import Qt, QTimer
import uuid
def merge_configs(general_options, specific_options):
    """
    Merges two dictionaries, with values from the `specific_options` dictionary
    overriding those in the `general_options` dictionary.

    Args:
        general_options: A dictionary containing general configuration options.
        specific_options: A dictionary containing specific configuration options.

    Returns:
        A dictionary containing the merged configuration options.
    """
    config = general_options.copy()
    config.update(specific_options)
    return config

class AnalysisWidgetBase(QWidget):
    """
    A base class for analysis widgets that provides common functionality such as 
    input validation, configuration merging, and job management.

    Attributes:
        job_manager: An instance of the JobManager responsible for managing job execution.
        timer: A QTimer instance that periodically checks the job manager queue.

    Methods:
        validate_inputs: Validates user inputs. Should be overridden by subclasses.
        run_analysis: Validates inputs, merges configurations, and runs the analysis.
        get_specific_options: Returns specific options for the analysis. Should be overridden by subclasses.
        run_specific_analysis: Runs the specific analysis. Should be overridden by subclasses.
        show_error_message: Displays an error message dialog.
        show_info_message: Displays an informational message dialog.
        on_job_finished: Handles the job finished event, displaying a message to the user.
        check_job_manager_queue: Periodically checks the job manager queue for updates.
    """
    def __init__(self, job_manager, general_options_getter=False):
        """
        Initializes the AnalysisWidgetBase with a job manager and an optional general options getter.

        Args:
            job_manager: The JobManager responsible for managing job execution.
            general_options_getter: An optional callable to retrieve general analysis options.
        """
        super().__init__()
        self.job_manager = job_manager
        self.job_manager.job_finished.connect(self.on_job_finished)
        self.setStyleSheet("""
            QToolBar {
                background-color: #333333;
                color: white;
                padding: 10px;
                font-size: 16px;
            }
            QPushButton {
                background-color: #555555;
                color: white;
                border: none;
                padding: 8px 16px;
                border-radius: 4px;
                margin: 0 5px;
                font-size: 16px;
            }
            QPushButton:hover {
                background-color: #666666;
            }
            QPushButton:pressed {
                background-color: #777777;
            }
        """)

        self.timer = QTimer()
        self.timer.timeout.connect(self.check_job_manager_queue)
        self.timer.start(1000)  # Check every second

    def validate_inputs(self):
        """
        Validates user inputs before running the analysis. 
        This method should be overridden by subclasses to provide specific validation logic.

        Returns:
            A list of error messages if validation fails, otherwise an empty list.
        """
        return []

    def run_analysis(self):
        """
        Validates inputs, merges general and specific configurations, and runs the analysis.
        If validation fails, an error message is displayed.
        """
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
        """
        Retrieves specific options for the analysis.
        This method should be overridden by subclasses to provide specific options.

        Returns:
            A dictionary containing specific options for the analysis.
        """
        return {}

    def run_specific_analysis(self, config):
        """
        Runs the specific analysis using the provided configuration.
        This method should be overridden by subclasses to implement specific analysis logic.

        Args:
            config: A dictionary containing the merged configuration options.
        """
        pass

    def show_error_message(self, errors):
        """
        Displays an error message dialog with a list of errors.

        Args:
            errors: A list of error messages to display.
        """
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Invalid input")
        msg.setInformativeText("\n".join(errors))
        msg.setWindowTitle("Error")
        msg.exec_()

    def show_info_message(self, message):
        """
        Displays an informational message dialog.

        Args:
            message: The message to display.
        """
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(message)
        msg.setWindowTitle("Info")
        msg.exec_()

    def on_job_finished(self, job_id, success, message):
        """
        Handles the job finished event, displaying an informational message to the user.

        Args:
            job_id: The ID of the finished job.
            success: A boolean indicating whether the job was successful.
            message: A message describing the outcome of the job.
        """
        status = "completed" if success else "failed"
        self.show_info_message(f"Job {job_id} {status}: {message}")

    def check_job_manager_queue(self):
        """
        Periodically checks the job manager queue for updates.
        This method can be overridden or extended if needed.
        """
        pass
