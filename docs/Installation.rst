Installing on Linux
=====================
Open a terminal, run the installation script by typing:

```
bash ./install.sh
```
Running the installation script will open the GUI.

To run it again, after exiting out of the GUI, type:

```
fast_ensemble/fast_ensemble-conda/bin/python -m fast_ensemble.gui.run_gui
```

For more information on usage, visit the Usage section.

Installing WSL
==============

To run FastEnsemble on Windows, make sure that WSL and Miniconda have been installed. Tutorials for installing are as follows. 

For the commands below, you need Windows 10 version 2004 and higher.

Open PowerShell, Windows Command Prompt, or Terminal in administrator mode. Do this by searching for one of the aforementioned in the search bar on the taskbar, and then right-clicking the option. Select “Run as administrator.”

Type the following command:

```
wsl --install
```

Then restart the machine.

Running FastEnsemble on Windows
===============================

Download the `install.sh` file from the GitHub repository.

Open the terminal and type `WSL` to enter the WSL environment.

Type `cd ~/` to navigate to the WSL file system.

Move the downloaded `install.sh` file to this directory by typing:

```
mv /location/of/file/install.sh ~/
```

Replace `/location/of/file/` with the path to where your `install.sh` file is located.

Run the installation script by typing:

```
bash ./install.sh
```

Running the installation script will open the GUI.

To run it again, after exiting out of the GUI, type:

```
fast_ensemble/fast_ensemble-conda/bin/python -m fast_ensemble.gui.run_gui
```

Note that you must be in the WSL file system (`~/`) to run this command successfully.

For more information on usage, visit the Usage section.
