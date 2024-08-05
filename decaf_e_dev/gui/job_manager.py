import multiprocessing as mp
from multiprocessing import Manager, Queue
from PyQt5.QtCore import pyqtSignal, QObject, Qt, pyqtSlot, QMetaObject
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QDialog, QTextEdit, QHBoxLayout
import threading
import uuid
import io
import sys
import pickle
import multiprocessing as mp
from multiprocessing import Manager, Queue
import uuid
import sys
import io

from PyQt5.QtCore import pyqtSignal, QObject
import threading
import multiprocessing as mp
from multiprocessing import Manager, Queue
import uuid
import sys
import io
import multiprocessing as mp
from multiprocessing import Manager, Queue
import uuid
import sys
import io


def job_wrapper(job_id, target, log_list, queue, *args):
    # Redirect stdout and stderr to capture logs
    sys.stdout = io.TextIOWrapper(io.BytesIO(), write_through=True)
    sys.stderr = sys.stdout
    try:
        target(*args)
        success = True
        message = "Job completed successfully."
    except Exception as e:
        success = False
        message = str(e)
    finally:
        log_list.append(sys.stdout.buffer.getvalue().decode())
        queue.put((job_id, success, message))
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

class JobManagerBackend:
    def __init__(self):
        self.jobs = {}
        self.manager = Manager()
        self.queue = self.manager.Queue()

    def run_job(self, target, args, job_name):
        # Ensure job_name is a string
        if not isinstance(job_name, str):
            raise ValueError("job_name must be a string")

        job_id = str(uuid.uuid4())
        if job_id in self.jobs:
            raise ValueError(f"Job ID {job_id} already exists.")

        # Check if args are pickleable
        for arg in args:
            try:
                pickle.dumps(arg)
            except pickle.PickleError as e:
                raise ValueError(f"Argument {arg} is not pickleable: {e}")
            except Exception as e:
                raise ValueError(f"Unexpected error while checking pickleability of {arg}: {e}")

            # Additional check for dictionary contents
            if isinstance(arg, dict):
                for key, value in arg.items():
                    try:
                        pickle.dumps(value)
                    except pickle.PickleError as e:
                        raise ValueError(f"Value for key '{key}' is not pickleable: {e}")
                    except Exception as e:
                        raise ValueError(f"Unexpected error while checking pickleability of value for key '{key}': {e}")

        log_list = self.manager.list()
        process = mp.Process(target=job_wrapper, args=(job_id, target, log_list, self.queue) + args)
        
        # Store only pickleable objects
        self.jobs[job_id] = {
            'name': job_name,
            'status': 'running',
            'log_list': log_list,
            'log': ''
        }
        process.start()
        return job_id

    def get_job_status(self, job_id):
        if job_id not in self.jobs:
            return None
        return self.jobs[job_id]['status']
    
    def get_job_name(self, job_id):
        if job_id not in self.jobs:
            return None
        return self.jobs[job_id]['name']

    def get_job_log(self, job_id):
        if job_id not in self.jobs:
            return None
        self.jobs[job_id]['log'] = ''.join(self.jobs[job_id]['log_list'])
        return self.jobs[job_id]['log']

    def update_job_status(self, job_id, status):
        if job_id in self.jobs:
            self.jobs[job_id]['status'] = status

from PyQt5.QtCore import pyqtSignal, QObject
import threading

class JobManager(QObject):
    job_finished = pyqtSignal(str, bool, str)  # job_id, success, message
    job_started = pyqtSignal(str)  # job_name

    def __init__(self):
        super().__init__()
        self.backend = JobManagerBackend()

        # Start a thread to monitor the queue
        self.monitor_thread = threading.Thread(target=self.monitor_queue, daemon=True)
        self.monitor_thread.start()

    def run_job(self, target, args, job_name):
        job_id = self.backend.run_job(target, args, job_name)
        self.job_started.emit(job_name)
        return job_id

    def get_job_status(self, job_id):
        return self.backend.get_job_status(job_id)

    def get_job_log(self, job_id):
        return self.backend.get_job_log(job_id)
    
    def get_job_name(self, job_id):
        return self.backend.get_job_name(job_id)

    def monitor_queue(self):
        while True:
            try:
                job_id, success, message = self.backend.queue.get()
                self.backend.update_job_status(job_id, "completed" if success else "failed")
                self.update_job_status(job_id, success, message)
            except EOFError:
                break

    def update_job_status(self, job_id, success, message):
        self.job_finished.emit(job_id, success, message)

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QDialog, QTextEdit, QHBoxLayout
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
        self.refresh_job_statuses()

    def refresh_job_statuses(self):
        self.list_widget.clear()
        for job_id, job_info in self.job_manager.backend.jobs.items():
            self.add_job_item(job_id, job_info['status'], "", job_info['name'])

    def update_job_status(self, job_id, success, message):
        status = "completed" if success else "failed"
        for index in range(self.list_widget.count()):
            item_widget = self.list_widget.itemWidget(self.list_widget.item(index))
            if item_widget and item_widget.job_id == job_id:
                item_widget.update_status(status, message)
                return
        self.add_job_item(job_id, status, message, self.job_manager.backend.jobs[job_id]['name'])

    def add_job_item(self, job_id, status, message, name):
        item_widget = JobItemWidget(job_id, status, message, self.job_manager, name)
        list_item = QListWidgetItem(self.list_widget)
        list_item.setSizeHint(item_widget.sizeHint())
        self.list_widget.addItem(list_item)
        self.list_widget.setItemWidget(list_item, item_widget)

class JobItemWidget(QWidget):
    def __init__(self, job_id, status, message, job_manager, name):
        super().__init__()
        self.job_id = job_id
        self.job_manager = job_manager
        self.name=name
        self.layout = QHBoxLayout(self)

        self.label = QLabel(f"Job {name} {status}: {message}")
        self.layout.addWidget(self.label)

        self.log_button = QPushButton("Log")
        self.log_button.clicked.connect(self.show_log)
        self.layout.addWidget(self.log_button)

    def update_status(self, status, message):
        self.label.setText(f"Job {self.name} {status}: {message}")

    def show_log(self):
        log_text = self.job_manager.get_job_log(self.job_id)
        log_dialog = QDialog(self)
        log_dialog.setWindowTitle(f"Job {self.name} Log")

        log_layout = QVBoxLayout(log_dialog)
        log_text_edit = QTextEdit()
        log_text_edit.setReadOnly(True)
        log_text_edit.setPlainText(log_text)
        log_layout.addWidget(log_text_edit)

        log_dialog.setLayout(log_layout)
        log_dialog.exec_()
