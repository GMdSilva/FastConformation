import multiprocessing as mp
from multiprocessing import Manager, Queue
from PyQt5.QtCore import pyqtSignal, QObject
import threading
import uuid

class JobManagerBackend:
    def __init__(self):
        self.jobs = {}
        self.queue = Queue()

    def run_job(self, target, args):
        job_id = str(uuid.uuid4())
        if job_id in self.jobs:
            raise ValueError(f"Job ID {job_id} already exists.")

        process = mp.Process(target=self._job_wrapper, args=(job_id, target) + args)
        self.jobs[job_id] = {
            'process': process,
            'status': 'running'
        }
        process.start()
        return job_id

    def _job_wrapper(self, job_id, target, *args):
        try:
            target(*args)
            success = True
            message = "Job completed successfully."
        except Exception as e:
            success = False
            message = str(e)
        finally:
            self.queue.put((job_id, success, message))

    def get_job_status(self, job_id):
        if job_id not in self.jobs:
            return None
        return self.jobs[job_id]['status']

class JobManager(QObject):
    job_finished = pyqtSignal(str, bool, str)  # job_id, success, message

    def __init__(self):
        super().__init__()
        self.backend = JobManagerBackend()

        # Start a thread to monitor the queue
        self.monitor_thread = threading.Thread(target=self.monitor_queue, daemon=True)
        self.monitor_thread.start()

    def run_job(self, target, args):
        return self.backend.run_job(target, args)

    def get_job_status(self, job_id):
        return self.backend.get_job_status(job_id)

    def monitor_queue(self):
        while True:
            job_id, success, message = self.backend.queue.get()
            self.update_job_status(job_id, success, message)

    def update_job_status(self, job_id, success, message):
        self.job_finished.emit(job_id, success, message)
        
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QListWidget, QListWidgetItem
from PyQt5.QtCore import Qt


class JobStatusPage(QWidget):
    def __init__(self, job_manager):
        super().__init__()
        self.job_manager = job_manager
        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignTop | Qt.AlignHCenter)

        self.title = QLabel("Job Status")
        font = self.title.font()
        font.setPointSize(20)
        self.title.setFont(font)
        self.title.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(self.title)

        self.list_widget = QListWidget(self)
        self.layout.addWidget(self.list_widget)

        self.job_manager.job_finished.connect(self.update_job_status)

    def update_job_status(self, job_id, success, message):
        status = "completed" if success else "failed"
        item = QListWidgetItem(f"Job {job_id} {status}: {message}")
        self.list_widget.addItem(item)
