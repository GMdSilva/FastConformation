import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QToolBar, QPushButton, QLineEdit, QComboBox, QLabel, QHBoxLayout, QCheckBox, QFileDialog
)
from PyQt5.QtCore import Qt, pyqtSignal, QThread
import uuid
from decaf_e_dev.predict_ensemble import run_ensemble_prediction

class JobStatusWidget(QWidget):
    analyze_signal = pyqtSignal(str)  # job_id

    def __init__(self, job_manager, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.job_manager = job_manager
        self.init_ui()

    def init_ui(self):
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.setWindowTitle("Job Status")

    def refresh_status(self):
        # Clear the layout
        for i in reversed(range(self.layout.count())): 
            widget = self.layout.itemAt(i).widget()
            if widget:
                widget.setParent(None)
        
        # Display all jobs
        for job in self.job_manager.get_jobs():
            job_layout = QHBoxLayout()
            job_label = QLabel(f"Job ID: {job['job_id']}, Type: {job['job_type']}, Status: {job['status']}")
            job_layout.addWidget(job_label)

            if job['status'] == 'Finished':
                analyze_button = QPushButton("Analyze Results")
                analyze_button.clicked.connect(lambda _, jid=job['job_id']: self.analyze_signal.emit(jid))
                job_layout.addWidget(analyze_button)
            elif job['status'] == 'Failed':
                error_label = QLabel(f"Error: {job['error']}")
                job_layout.addWidget(error_label)

            self.layout.addLayout(job_layout)

    def update_status(self, job_id, message):
        self.job_manager.update_status(job_id, message)
        self.refresh_status()

    def job_finished(self, job_id):
        self.refresh_status()
class JobManager:
    def __init__(self):
        self.jobs = []

    def add_job(self, job_id, job_type, config):
        self.jobs.append({
            'job_id': job_id,
            'job_type': job_type,
            'config': config,
            'status': 'Pending',
            'error': None
        })

    def update_status(self, job_id, status, error=None):
        for job in self.jobs:
            if job['job_id'] == job_id:
                job['status'] = status
                job['error'] = error
                break

    def get_jobs(self):
        return self.jobs
class WorkerThread(QThread):
    progress = pyqtSignal(str, str)  # job_id, message
    finished = pyqtSignal(str)  # job_id

    def __init__(self, job_id, job_type, config):
        super().__init__()
        self.job_id = job_id
        self.job_type = job_type
        self.config = config

    def run(self):
        try:
            self.progress.emit(self.job_id, "Starting job...")
            run_ensemble_prediction(self.config)
            self.progress.emit(self.job_id, "Job finished successfully.")
        except Exception as e:
            self.progress.emit(self.job_id, f"Job failed: {str(e)}")
        self.finished.emit(self.job_id)
