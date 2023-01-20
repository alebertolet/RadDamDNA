#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 1/20/23 9:23 AM

@author: alejandrobertolet
"""
import sys, io
from PyQt5 import QtCore, QtGui, QtWidgets
from RadDamDNA.bioStage.running import Simulator

class SimulatorThread(QtCore.QThread):
    def __init__(self, time_options, diffusion_model, dsb_model, ssb_model, bd_model, nucleus_max_radius, basepath, max_dose, version, n_runs):
        super().__init__()
        self.time_options = time_options
        self.diffusion_model = diffusion_model
        self.dsb_model = dsb_model
        self.ssb_model = ssb_model
        self.bd_model = bd_model
        self.nucleus_max_radius = nucleus_max_radius
        self.basepath = basepath
        self.max_dose = max_dose
        self.version = version
        self.n_runs = n_runs

    def run(self):
        sim = Simulator(self.time_options, self.diffusion_model, self.dsb_model, self.ssb_model, self.bd_model, self.nucleus_max_radius)
        sim.ReadDamage(self.basepath, self.max_dose, self.version)
        sim.Run(self.n_runs, rereadDamageForNewRuns=False, basepath=self.basepath, maxDose=self.max_dose, version=self.version)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        # Create widgets for inputting parameters
        self.initial_time_label = QtWidgets.QLabel("Initial Time (h):")
        self.initial_time_input = QtWidgets.QLineEdit()
        self.final_time_label = QtWidgets.QLabel("Final Time (h):")
        self.final_time_input = QtWidgets.QLineEdit()
        self.time_points_label = QtWidgets.QLabel("Time Points:")
        self.time_points_input = QtWidgets.QLineEdit()
        self.nucleus_max_radius_label = QtWidgets.QLabel("Nucleus Max Radius:")
        self.nucleus_max_radius_input = QtWidgets.QLineEdit()
        self.diffusion_model_label = QtWidgets.QLabel("Diffusion Model:")
        self.diffusion_model_input = QtWidgets.QComboBox()
        self.diffusion_model_input.addItems(["free", "none"])
        self.dsb_model_label = QtWidgets.QLabel("DSB repair Model:")
        self.dsb_model_input = QtWidgets.QComboBox()
        self.dsb_model_input.addItems(["standard", "none"])
        self.ssb_model_label = QtWidgets.QLabel("SSB repair Model:")
        self.ssb_model_input = QtWidgets.QComboBox()
        self.ssb_model_input.addItems(["standard", "none"])
        self.bd_model_label = QtWidgets.QLabel("BD repair Model:")
        self.bd_model_input = QtWidgets.QComboBox()
        self.bd_model_input.addItems(["standard", "none"])
        self.n_runs_label = QtWidgets.QLabel("Number of Runs:")
        self.n_runs_input = QtWidgets.QLineEdit()
        self.basepath_label = QtWidgets.QLabel("Basepath:")
        self.basepath_input = QtWidgets.QLineEdit()
        self.basepath_browse_button = QtWidgets.QPushButton("Browse")
        self.basepath_browse_button.clicked.connect(self.browse_basepath)
        self.max_dose_label = QtWidgets.QLabel("Max Dose:")
        self.max_dose_input = QtWidgets.QLineEdit()
        self.version_label = QtWidgets.QLabel("Version:")
        self.version_input = QtWidgets.QLineEdit()
        self.run_button = QtWidgets.QPushButton("Run")

        # Default values
        self.initial_time_input.setText("0")
        self.final_time_input.setText("25")
        self.time_points_input.setText("100")
        self.nucleus_max_radius_input.setText("4.65")
        self.n_runs_input.setText("10")
        self.basepath_input.setText("/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/xray/sims/250keV.txt/")
        self.max_dose_input.setText("0.25")
        self.version_input.setText("1.0")

        # Create layout and add widgets
        self.layout = QtWidgets.QFormLayout()
        # Create layout and add widgets
        self.layout.addRow(self.initial_time_label, self.initial_time_input)
        self.layout.addRow(self.final_time_label, self.final_time_input)
        self.layout.addRow(self.time_points_label, self.time_points_input)
        self.layout.addRow(self.nucleus_max_radius_label, self.nucleus_max_radius_input)
        self.layout.addRow(self.diffusion_model_label, self.diffusion_model_input)
        self.layout.addRow(self.dsb_model_label, self.dsb_model_input)
        self.layout.addRow(self.ssb_model_label, self.ssb_model_input)
        self.layout.addRow(self.bd_model_label, self.bd_model_input)
        self.layout.addRow(self.n_runs_label, self.n_runs_input)
        self.basepath_layout = QtWidgets.QHBoxLayout()
        self.basepath_layout.addWidget(self.basepath_input)
        self.basepath_layout.addWidget(self.basepath_browse_button)
        self.layout.addRow(self.basepath_label, self.basepath_layout)
        self.layout.addRow(self.max_dose_label, self.max_dose_input)
        self.layout.addRow(self.version_label, self.version_input)
        self.layout.addRow(self.run_button)

        # Create central widget and set layout
        self.central_widget = QtWidgets.QWidget()
        self.central_widget.setLayout(self.layout)
        self.setCentralWidget(self.central_widget)

        # Connect run button to run function
        self.run_button.clicked.connect(self.run)

    def browse_basepath(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.ReadOnly
        directory = QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory", options=options)
        if directory:
            self.basepath_input.setText(directory)

    def run(self):
        # Get input values
        initial_time = int(self.initial_time_input.text()) * 3600
        final_time = int(self.final_time_input.text()) * 3600
        time_points = int(self.time_points_input.text())
        time_options = [initial_time, final_time, time_points]
        nucleus_max_radius = float(self.nucleus_max_radius_input.text())
        diffusion_model = self.diffusion_model_input.currentText()
        dsb_model = self.dsb_model_input.currentText()
        ssb_model = self.ssb_model_input.currentText()
        bd_model = self.bd_model_input.currentText()
        n_runs = int(self.n_runs_input.text())
        basepath = self.basepath_input.text()
        max_dose = float(self.max_dose_input.text())
        version = self.version_input.text()

        # Create instance of Simulator and set parameters
        self.simulator_thread = SimulatorThread(time_options, diffusion_model, dsb_model, ssb_model, bd_model,
                                                nucleus_max_radius, basepath, max_dose, version, n_runs)
        self.simulator_thread.start()

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())