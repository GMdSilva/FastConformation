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

from queue import Empty

import io
import sys
import uuid
import pickle
import multiprocessing as mp
from multiprocessing import Manager
from queue import Empty
from PyQt5.QtCore import pyqtSignal, QObject, Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QDialog, QTextEdit, QHBoxLayout

def job_wrapper(job_id, target, log_list, queue, *args):
    """
    Wraps the execution of a job to capture stdout and stderr output into a log list, 
    and handles exceptions and job completion status.

    Args:
        job_id: A unique identifier for the job.
        target: The function to be executed as the job.
        log_list: A shared list to capture log messages.
        queue: A multiprocessing queue to report the job status.
        *args: Additional arguments to pass to the target function.
    """
    class StreamToLogger(io.StringIO):
        def __init__(self, log_list):
            super().__init__()
            self.log_list = log_list

        def write(self, message):
            super().write(message)
            self.log_list.append(message)

    stdout_logger = StreamToLogger(log_list)
    stderr_logger = StreamToLogger(log_list)
    sys.stdout = stdout_logger
    sys.stderr = stderr_logger

    try:
        target(*args)
        success = True
        message = "Job completed successfully."
    except KeyboardInterrupt:
        success = False
        message = "Job was interrupted by the user."
        raise  # Re-raise the exception to allow the program to terminate
    except Exception as e:
        success = False
        message = str(e)
    finally:
        queue.put((job_id, success, message))
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

class JobManagerBackend:
    """
    Manages the backend operations of job execution, including job tracking, status updates, 
    and log management.

    Attributes:
        jobs: A dictionary to store information about the jobs.
        manager: A multiprocessing.Manager instance to manage shared resources.
        queue: A multiprocessing.Queue for communication between processes.
    """
    def __init__(self):
        self.jobs = {}
        self.manager = Manager()
        self.queue = self.manager.Queue()

    def run_job(self, target, args, job_name):
        """
        Runs a job in a separate process, ensuring that the job name is unique and 
        that the arguments are pickleable.

        Args:
            target: The function to execute as the job.
            args: A tuple of arguments to pass to the target function.
            job_name: A string name for the job.

        Returns:
            job_id: A unique identifier for the job.

        Raises:
            ValueError: If job_name is not a string, or if arguments are not pickleable.
        """
        if not isinstance(job_name, str):
            raise ValueError("job_name must be a string")

        job_id = str(uuid.uuid4())
        if job_id in self.jobs:
            raise ValueError(f"Job ID {job_id} already exists.")

        for arg in args:
            try:
                pickle.dumps(arg)
            except pickle.PickleError as e:
                raise ValueError(f"Argument {arg} is not pickleable: {e}")
            except Exception as e:
                raise ValueError(f"Unexpected error while checking pickleability of {arg}: {e}")

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
        
        self.jobs[job_id] = {
            'name': job_name,
            'status': 'running',
            'log_list': log_list,
            'log': ''
        }
        process.start()
        return job_id

    def get_job_status(self, job_id):
        """
        Retrieves the current status of the specified job.

        Args:
            job_id: The unique identifier of the job.

        Returns:
            The status of the job as a string, or None if the job_id is not found.
        """
        if job_id not in self.jobs:
            return None
        return self.jobs[job_id]['status']
    
    def get_job_name(self, job_id):
        """
        Retrieves the name of the specified job.

        Args:
            job_id: The unique identifier of the job.

        Returns:
            The name of the job as a string, or None if the job_id is not found.
        """
        if job_id not in self.jobs:
            return None
        return self.jobs[job_id]['name']

    def get_job_log(self, job_id):
        """
        Retrieves the log of the specified job.

        Args:
            job_id: The unique identifier of the job.

        Returns:
            The log of the job as a string, or None if the job_id is not found.
        """
        if job_id not in self.jobs:
            return None
        self.jobs[job_id]['log'] = ''.join(self.jobs[job_id]['log_list'])
        return self.jobs[job_id]['log']

    def update_job_status(self, job_id, status):
        """
        Updates the status of the specified job.

        Args:
            job_id: The unique identifier of the job.
            status: The new status of the job as a string.
        """
        if job_id in self.jobs:
            self.jobs[job_id]['status'] = status

class JobManager(QObject):
    """
    Manages job execution, monitoring, and status updates within a Qt application. 
    It handles starting jobs, updating statuses, and emitting signals.

    Attributes:
        job_finished: A PyQt signal emitted when a job finishes, with job_id, success, and message.
        job_started: A PyQt signal emitted when a job starts, with the job name.
    """
    job_finished = pyqtSignal(str, bool, str)  # job_id, success, message
    job_started = pyqtSignal(str)  # job_name

    def __init__(self):
        super().__init__()
        self.backend = JobManagerBackend()
        self._stop_event = threading.Event()
        self.monitor_thread = threading.Thread(target=self.monitor_queue, daemon=True)
        self.monitor_thread.start()

    def run_job(self, target, args, job_name):
        """
        Runs a job using the backend and emits a signal when the job starts.

        Args:
            target: The function to execute as the job.
            args: A tuple of arguments to pass to the target function.
            job_name: A string name for the job.

        Returns:
            job_id: The unique identifier for the job.
        """
        job_id = self.backend.run_job(target, args, job_name)
        self.job_started.emit(job_name)
        return job_id

    def get_job_status(self, job_id):
        """
        Retrieves the current status of the specified job.

        Args:
            job_id: The unique identifier of the job.

        Returns:
            The status of the job as a string.
        """
        return self.backend.get_job_status(job_id)

    def get_job_log(self, job_id):
        """
        Retrieves the log of the specified job.

        Args:
            job_id: The unique identifier of the job.

        Returns:
            The log of the job as a string.
        """
        return self.backend.get_job_log(job_id)
    
    def get_job_name(self, job_id):
        """
        Retrieves the name of the specified job.

        Args:
            job_id: The unique identifier of the job.

        Returns:
            The name of the job as a string.
        """
        return self.backend.get_job_name(job_id)

    def monitor_queue(self):
        """
        Monitors the queue for job status updates and emits the job_finished signal 
        when a job is completed or fails.
        """
        while not self._stop_event.is_set():
            try:
                job_id, success, message = self.backend.queue.get(timeout=1)  # Add timeout to allow periodic check
                self.backend.update_job_status(job_id, "completed" if success else "failed")
                self.update_job_status(job_id, success, message)
            except Empty:
                continue
            
    def stop(self):
        """
        Stops the monitoring thread for job status updates.
        """
        self._stop_event.set()
        self.monitor_thread.join()

    def update_job_status(self, job_id, success, message):
        """
        Emits the job_finished signal with the job's status.

        Args:
            job_id: The unique identifier of the job.
            success: A boolean indicating if the job was successful.
            message: A string message about the job's outcome.
        """
        self.job_finished.emit(job_id, success, message)

class JobStatusPage(QWidget):
    """
    Displays the status of all jobs in the JobManager. It allows users to view 
    the current status and logs of jobs.

    Methods:
        refresh_job_statuses: Refreshes the list of job statuses.
        update_job_status: Updates the displayed status of a specific job.
        add_job_item: Adds a new job item to the list.
    """
    def __init__(self, job_manager):
        """
        Initialize the JobStatusPage with a JobManager instance.

        Args:
            job_manager: The JobManager responsible for managing job execution and status.
        """
        super().__init__()
        self.job_manager = job_manager
        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignTop | Qt.AlignHCenter)

        self.title = QLabel("Job Status")
        font = self.title.font()
        font.setPointSize(16)
        self.title.setFont(font)
        self.title.setAlignment(Qt.AlignCenter)
        self.layout.addWidget(self.title)

        self.list_widget = QListWidget(self)
        self.layout.addWidget(self.list_widget)
        self.job_manager.job_finished.connect(self.update_job_status)
        self.refresh_job_statuses()
        
    def refresh_job_statuses(self):
        """
        Clears and repopulates the list widget with the current job statuses.
        """
        self.list_widget.clear()
        for job_id, job_info in self.job_manager.backend.jobs.items():
            self.add_job_item(job_id, job_info['status'], "", job_info['name'])

    def update_job_status(self, job_id, success, message):
        """
        Updates the status of a specific job in the list widget.

        Args:
            job_id: The unique identifier of the job.
            success: A boolean indicating if the job was successful.
            message: A string message about the job's outcome.
        """
        status = "completed" if success else "failed"
        for index in range(self.list_widget.count()):
            item_widget = self.list_widget.itemWidget(self.list_widget.item(index))
            if item_widget and item_widget.job_id == job_id:
                item_widget.update_status(status, message)
                return
        self.add_job_item(job_id, status, message, self.job_manager.backend.jobs[job_id]['name'])

    def add_job_item(self, job_id, status, message, name):
        """
        Adds a new job item to the list widget.

        Args:
            job_id: The unique identifier of the job.
            status: The current status of the job as a string.
            message: A string message about the job's outcome.
            name: The name of the job.
        """
        item_widget = JobItemWidget(job_id, status, message, self.job_manager, name)
        list_item = QListWidgetItem(self.list_widget)
        list_item.setSizeHint(item_widget.sizeHint())
        self.list_widget.addItem(list_item)
        self.list_widget.setItemWidget(list_item, item_widget)

class JobItemWidget(QWidget):
    """
    Represents a single job item in the JobStatusPage, displaying the job's name, 
    status, and a button to view the job's log.

    Methods:
        update_status: Updates the displayed status and message of the job.
        show_log: Displays the log of the job in a new dialog.
    """
    def __init__(self, job_id, status, message, job_manager, name):
        """
        Initialize the JobItemWidget with job details.

        Args:
            job_id: The unique identifier of the job.
            status: The current status of the job as a string.
            message: A string message about the job's outcome.
            job_manager: The JobManager responsible for managing job execution and status.
            name: The name of the job.
        """
        super().__init__()
        self.job_id = job_id
        self.job_manager = job_manager
        self.name = name
        self.layout = QHBoxLayout(self)

        self.label = QLabel(f"Job {name} {status}: {message}")
        self.layout.addWidget(self.label)

        self.log_button = QPushButton("Show Log")
        self.log_button.setMinimumHeight(50)
        self.log_button.clicked.connect(self.show_log)
        self.layout.addWidget(self.log_button)

    def update_status(self, status, message):
        """
        Updates the displayed status and message of the job.

        Args:
            status: The new status of the job as a string.
            message: A string message about the job's outcome.
        """
        self.label.setText(f"Job {self.name} {status}: {message}")

    def show_log(self):
        """
        Displays the log of the job in a new dialog window.
        """
        log_text = self.job_manager.get_job_log(self.job_id)
        log_dialog = QDialog(self)
        log_dialog.setWindowTitle(f"Job {self.name} Log")
        log_dialog.setMinimumWidth(500)

        log_layout = QVBoxLayout(log_dialog)
        log_text_edit = QTextEdit()
        log_text_edit.setReadOnly(True)
        log_text_edit.setPlainText(log_text)
        
        log_layout.addWidget(log_text_edit)

        log_dialog.setLayout(log_layout)
        log_dialog.exec_()
