from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QFormLayout, QLineEdit,
    QPushButton, QMessageBox, QSpinBox, QDoubleSpinBox, QTabWidget,
    QRadioButton, QButtonGroup, QHBoxLayout, QFileDialog, QCheckBox, QLabel,
    QComboBox, QScrollArea
)
import sys
import json5
import json
import subprocess
from PySide6.QtGui import QIcon
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import QDialog, QVBoxLayout, QProgressBar, QPushButton, QLabel, QPlainTextEdit
from PySide6.QtGui import QFontDatabase

import subprocess
import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
matplotlib.rcParams["font.family"] = "Serif"
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT

from pyvistaqt import QtInteractor


class SolverThread(QThread):
    output = Signal(str)  # New signal for live output
    finished = Signal(int, str, str)

    def __init__(self, cmd):
        super().__init__()
        self.cmd = cmd
        self.proc = None

    def run(self):
        self.proc = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        output_lines = []
        for line in self.proc.stdout:
            self.output.emit(line)
            output_lines.append(line)
        self.proc.stdout.close()
        self.proc.wait()
        self.finished.emit(self.proc.returncode, ''.join(output_lines), '')  # stderr is merged into stdout

    def terminate_solver(self):
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()

class ProgressDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Running Solver...")
        layout = QVBoxLayout(self)
        self.label = QLabel("Solver is running...")
        layout.addWidget(self.label)
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)  # Indeterminate
        layout.addWidget(self.progress)
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        # Set fixed-width font
        self.console.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        layout.addWidget(self.console)
        self.cancel_btn = QPushButton("Cancel")
        layout.addWidget(self.cancel_btn)

class HaloConfigGUI(QWidget):
    def __init__(self):
        super().__init__()

        # self.setStyleSheet("QLineEdit { min-width: 0; width: 100%; }")
        # self.setStyleSheet("QLineEdit { min-width: 0; }")
        self.setWindowTitle("Halo Solver")

        # 1. Create the main vertical layout for the whole window
        self.outer_layout = QVBoxLayout(self)
        self.setLayout(self.outer_layout)

        # 2. Create the main horizontal layout (tabs + plot)
        self.main_layout = QHBoxLayout()
        self.outer_layout.addLayout(self.main_layout, stretch=1)

        # 3. Add tabs and plot panel to self.main_layout as before
        self.tabs = QTabWidget()
        self.main_layout.addWidget(self.tabs, stretch=2)

        # Helper to wrap a widget in a scroll area
        def make_scrollable(widget):
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(widget)
            return scroll

        # --- Epoch Tab ---
        self.epoch_tab = QWidget()
        self.epoch_form = QFormLayout(self.epoch_tab)
        self.radio_group = QButtonGroup(self)
        self.radio_calendar = QRadioButton("Calendar Date")
        self.radio_et = QRadioButton("Ephemeris Time")
        self.radio_calendar.setChecked(True)
        self.radio_group.addButton(self.radio_calendar)
        self.radio_group.addButton(self.radio_et)
        radio_layout = QHBoxLayout()
        radio_layout.addWidget(self.radio_calendar)
        radio_layout.addWidget(self.radio_et)
        self.epoch_form.addRow("Epoch Mode (TDB):", radio_layout)
        self.year = QSpinBox(); self.year.setRange(1900, 2100); self.year.setValue(2025)
        self.month = QSpinBox(); self.month.setRange(1, 12)
        self.day = QSpinBox(); self.day.setRange(1, 31)
        self.hour = QSpinBox(); self.hour.setRange(0, 23)
        self.minute = QSpinBox(); self.minute.setRange(0, 59)
        self.sec = QDoubleSpinBox(); self.sec.setRange(0, 59.999); self.sec.setDecimals(3)
        self.epoch_form.addRow("Year:", self.year)
        self.epoch_form.addRow("Month:", self.month)
        self.epoch_form.addRow("Day:", self.day)
        self.epoch_form.addRow("Hour:", self.hour)
        self.epoch_form.addRow("Minute:", self.minute)
        self.epoch_form.addRow("Second:", self.sec)
        self.et_ref = QDoubleSpinBox(); self.et_ref.setDecimals(6)
        self.et_ref.setRange(-1e12, 1e12)
        self.et_ref.setValue(0.0)
        self.epoch_form.addRow("Ephemeris Time (sec):", self.et_ref)
        self.tabs.addTab(make_scrollable(self.epoch_tab), "Epoch")

        # --- Mission Tab ---
        self.mission_tab = QWidget()
        self.mission_form = QFormLayout(self.mission_tab)
        self.N_or_S = QLineEdit("N")
        self.L1_or_L2 = QLineEdit("L1")
        self.period = QDoubleSpinBox(); self.period.setDecimals(10); self.period.setRange(0, 100)
        self.n_revs = QSpinBox(); self.n_revs.setRange(0, 10000)
        self.mission_form.addRow("N_or_S:", self.N_or_S)
        self.mission_form.addRow("L1_or_L2:", self.L1_or_L2)
        self.mission_form.addRow("period:", self.period)
        self.mission_form.addRow("n_revs:", self.n_revs)

        self.patch_point_file = QLineEdit("data/L2_halos.json")
        self.patch_point_file_is_periapsis = QCheckBox()
        self.initial_guess_from_file = QLineEdit()
        self.mission_form.addRow("patch_point_file:", self.patch_point_file)
        self.mission_form.addRow("patch_point_file_is_periapsis:", self.patch_point_file_is_periapsis)
        self.mission_form.addRow("initial_guess_from_file:", self.initial_guess_from_file)

        self.tabs.addTab(make_scrollable(self.mission_tab), "Mission")

        # --- Optimizer Tab ---
        self.sa_tab = QWidget()
        self.sa_form = QFormLayout(self.sa_tab)

        self.sa_lb_0 = QDoubleSpinBox(); self.sa_lb_0.setDecimals(6); self.sa_lb_0.setRange(-1e12, 1e12)
        self.sa_lb_1 = QDoubleSpinBox(); self.sa_lb_1.setDecimals(6); self.sa_lb_1.setRange(-1e12, 1e12)
        self.sa_ub_0 = QDoubleSpinBox(); self.sa_ub_0.setDecimals(6); self.sa_ub_0.setRange(-1e12, 1e12)
        self.sa_ub_1 = QDoubleSpinBox(); self.sa_ub_1.setDecimals(6); self.sa_ub_1.setRange(-1e12, 1e12)
        self.vm_0 = QDoubleSpinBox(); self.vm_0.setDecimals(6); self.vm_0.setRange(-1e12, 1e12)
        self.vm_1 = QDoubleSpinBox(); self.vm_1.setDecimals(6); self.vm_1.setRange(-1e12, 1e12)
        self.sa_ns = QSpinBox(); self.sa_ns.setRange(0, 1000)
        self.sa_nt = QSpinBox(); self.sa_nt.setRange(0, 1000)
        self.sa_neps = QSpinBox(); self.sa_neps.setRange(0, 1000)
        self.sa_form.addRow("sa_lb[0]:", self.sa_lb_0)
        self.sa_form.addRow("sa_lb[1]:", self.sa_lb_1)
        self.sa_form.addRow("sa_ub[0]:", self.sa_ub_0)
        self.sa_form.addRow("sa_ub[1]:", self.sa_ub_1)
        self.sa_form.addRow("vm[0]:", self.vm_0)
        self.sa_form.addRow("vm[1]:", self.vm_1)
        self.sa_form.addRow("sa_ns:", self.sa_ns)
        self.sa_form.addRow("sa_nt:", self.sa_nt)
        self.sa_form.addRow("sa_neps:", self.sa_neps)

        self.tabs.addTab(make_scrollable(self.sa_tab), "Optimizer")

        # --- Solver Tab ---
        self.solver_tab = QWidget()
        self.solver_form = QFormLayout(self.solver_tab)

        # self.solver_form.addRow("solver_mode:", self.solver_mode)
        # Solver mode menu
        self.solver_mode_menu = QComboBox()
        self.solver_mode_options = [
            ("Dense solver (LAPACK)", 1),
            ("Sparse solver (LSQR)", 2),
            ("Sparse solver (LUSOL)", 3),
            ("Sparse solver (LMSR)", 4),
            ("Custom solver (QR_MUMPS)", 5)
        ]
        for label, value in self.solver_mode_options:
            self.solver_mode_menu.addItem(label, value)
        self.solver_form.addRow("solver_mode:", self.solver_mode_menu)
        self.nlesolver_tol = QLineEdit("1.0e-6")
        self.solver_form.addRow("nlesolver_tol:", self.nlesolver_tol)

        self.fix_initial_time = QCheckBox(); self.fix_initial_time.setChecked(False)
        self.fix_initial_r = QCheckBox(); self.fix_initial_r.setChecked(False)
        self.fix_ry_at_end_of_rev = QSpinBox(); self.fix_ry_at_end_of_rev.setRange(-1000, 1000)
        self.fix_final_ry_and_vx = QCheckBox(); self.fix_final_ry_and_vx.setChecked(False)
        self.constrain_initial_rdot = QCheckBox(); self.constrain_initial_rdot.setChecked(False)
        self.fscale_rdot = QDoubleSpinBox(); self.fscale_rdot.setDecimals(6); self.fscale_rdot.setRange(-1e12, 1e12)
        self.solver_mode = QSpinBox(); self.solver_mode.setRange(1, 10)
        self.solver_form.addRow("fix_initial_time:", self.fix_initial_time)
        self.solver_form.addRow("fix_initial_r:", self.fix_initial_r)
        self.solver_form.addRow("fix_ry_at_end_of_rev:", self.fix_ry_at_end_of_rev)
        self.solver_form.addRow("fix_final_ry_and_vx:", self.fix_final_ry_and_vx)
        self.solver_form.addRow("constrain_initial_rdot:", self.constrain_initial_rdot)
        self.solver_form.addRow("fscale_rdot:", self.fscale_rdot)

        self.tabs.addTab(make_scrollable(self.solver_tab), "Solver")


        # --- Integrator Tab ---
        self.integrator_tab = QWidget()
        self.integrator_form = QFormLayout(self.integrator_tab)
        self.rtol = QLineEdit("1.0e-12")
        self.atol = QLineEdit("1.0e-12")
        self.integrator_form.addRow("rtol:", self.rtol)
        self.integrator_form.addRow("atol:", self.atol)
        self.tabs.addTab(make_scrollable(self.integrator_tab), "Integrator")

        # --- Force Model Tab ---
        self.force_tab = QWidget()
        self.force_form = QFormLayout(self.force_tab)
        self.use_splined_ephemeris = QCheckBox(); self.use_splined_ephemeris.setChecked(True)
        self.dt_spline_sec = QSpinBox(); self.dt_spline_sec.setRange(1, 1_000_000)
        self.dt_spline_sec.setValue(3600)
        self.ephemeris_file = QLineEdit("data/eph/JPLEPH.421")
        self.gravfile = QLineEdit("data/grav/gggrx_0020pm_sha.tab")
        self.grav_n = QSpinBox(); self.grav_n.setRange(0, 1000)
        self.grav_m = QSpinBox(); self.grav_m.setRange(0, 1000)
        # self.grav_frame = QSpinBox(); self.grav_frame.setRange(0, 10)
        self.grav_frame = QComboBox()
        self.grav_frame.addItem("iau_moon", 1)
        self.grav_frame.addItem("moon_pa", 2)
        self.moon_pa_file = QLineEdit("data/moon_pa_2000_2100.csv")

        self.pointmass_central_body = QCheckBox()
        self.include_pointmass_earth = QCheckBox(); self.include_pointmass_earth.setChecked(True)
        self.include_pointmass_sun = QCheckBox(); self.include_pointmass_sun.setChecked(True)
        self.include_pointmass_jupiter = QCheckBox()
        self.use_battin_gravity = QCheckBox()
        self.force_form.addRow("use_splined_ephemeris:", self.use_splined_ephemeris)
        self.force_form.addRow("dt_spline_sec:", self.dt_spline_sec)
        self.force_form.addRow("ephemeris_file:", self.ephemeris_file)
        self.force_form.addRow("gravfile:", self.gravfile)
        self.force_form.addRow("grav_n:", self.grav_n)
        self.force_form.addRow("grav_m:", self.grav_m)
        self.force_form.addRow("grav_frame:", self.grav_frame)
        self.force_form.addRow("moon_pa_file:", self.moon_pa_file)
        self.force_form.addRow("pointmass_central_body:", self.pointmass_central_body)
        self.force_form.addRow("include_pointmass_earth:", self.include_pointmass_earth)
        self.force_form.addRow("include_pointmass_sun:", self.include_pointmass_sun)
        self.force_form.addRow("include_pointmass_jupiter:", self.include_pointmass_jupiter)
        self.force_form.addRow("use_battin_gravity:", self.use_battin_gravity)
        self.tabs.addTab(make_scrollable(self.force_tab), "Force Model")

        # --- Output Tab ---
        self.output_tab = QWidget()
        self.output_form = QFormLayout(self.output_tab)
        self.generate_plots = QCheckBox()
        self.generate_trajectory_files = QCheckBox(); self.generate_trajectory_files.setChecked(True)
        self.generate_kernel = QCheckBox(); self.generate_kernel.setChecked(True)
        self.object_id = QSpinBox(); self.object_id.setRange(-1_000_000, 1_000_000)
        self.object_id.setValue(-50000)
        self.object_name = QLineEdit("HALO")
        self.leapseconds_file = QLineEdit("kernel/naif0012.tls")
        self.mkspk_path = QLineEdit("kernel/mkspk")
        self.polynom_degree = QSpinBox(); self.polynom_degree.setRange(0, 20)
        self.polynom_degree.setValue(9)
        self.output_spk_type = QSpinBox(); self.output_spk_type.setRange(0, 20)
        self.output_spk_type.setValue(9)
        self.segment_id = QLineEdit("SPK_STATES_09")
        self.generate_guess_and_solution_files = QCheckBox(); self.generate_guess_and_solution_files.setChecked(True)
        self.generate_defect_file = QCheckBox(); self.generate_defect_file.setChecked(True)
        self.generate_json_trajectory_file = QCheckBox(); self.generate_json_trajectory_file.setChecked(True)
        self.generate_eclipse_files = QCheckBox(); self.generate_eclipse_files.setChecked(True)
        self.r_eclipse_bubble = QDoubleSpinBox(); self.r_eclipse_bubble.setDecimals(3); self.r_eclipse_bubble.setRange(0, 1e6)
        self.r_eclipse_bubble.setValue(1.0)
        self.eclipse_dt_step = QDoubleSpinBox(); self.eclipse_dt_step.setDecimals(3); self.eclipse_dt_step.setRange(0, 1e7)
        self.eclipse_dt_step.setValue(3600.0)
        self.eclipse_filetype = QSpinBox(); self.eclipse_filetype.setRange(0, 10)
        self.eclipse_filetype.setValue(2)
        self.run_pyvista_script = QCheckBox(); self.run_pyvista_script.setChecked(True)
        self.generate_rp_ra_file = QCheckBox(); self.generate_rp_ra_file.setChecked(True)
        self.output_form.addRow("generate_plots:", self.generate_plots)
        self.output_form.addRow("generate_trajectory_files:", self.generate_trajectory_files)
        self.output_form.addRow("generate_kernel:", self.generate_kernel)
        self.output_form.addRow("object_id:", self.object_id)
        self.output_form.addRow("object_name:", self.object_name)
        self.output_form.addRow("leapseconds_file:", self.leapseconds_file)
        self.output_form.addRow("mkspk_path:", self.mkspk_path)
        self.output_form.addRow("polynom_degree:", self.polynom_degree)
        self.output_form.addRow("output_spk_type:", self.output_spk_type)
        self.output_form.addRow("segment_id:", self.segment_id)
        self.output_form.addRow("generate_guess_and_solution_files:", self.generate_guess_and_solution_files)
        self.output_form.addRow("generate_defect_file:", self.generate_defect_file)
        self.output_form.addRow("generate_json_trajectory_file:", self.generate_json_trajectory_file)
        self.output_form.addRow("generate_eclipse_files:", self.generate_eclipse_files)
        self.output_form.addRow("r_eclipse_bubble:", self.r_eclipse_bubble)
        self.output_form.addRow("eclipse_dt_step:", self.eclipse_dt_step)
        self.output_form.addRow("eclipse_filetype:", self.eclipse_filetype)
        self.output_form.addRow("run_pyvista_script:", self.run_pyvista_script)
        self.output_form.addRow("generate_rp_ra_file:", self.generate_rp_ra_file)
        self.tabs.addTab(make_scrollable(self.output_tab), "Output")

        # Matplotlib plot panel (right)
        # self.plot_panel = QWidget()
        # self.plot_layout = QVBoxLayout(self.plot_panel)
        # self.figure = Figure(figsize=(5, 4))
        # self.canvas = FigureCanvas(self.figure)
        # self.plot_layout.addWidget(self.canvas)
        # self.main_layout.addWidget(self.plot_panel, stretch=3)

        # PyVista plot panel (right)
        #self.plotter = QtInteractor(self, shape=(1,3))  # to use the subplots
        self.plotter = QtInteractor(self, shape=(1,1))
        self.main_layout.addWidget(self.plotter, stretch=3)

        # self.toolbar = NavigationToolbar2QT(self.canvas, self.plot_panel)
        # self.plot_layout.addWidget(self.toolbar)
        self.ax = None

        # Example: plot something initially
        self.plot_example()

        # --- File name field and buttons ---
        filename_layout = QHBoxLayout()
        self.filename = QLineEdit("config.json")
        filename_layout.addWidget(QLabel("Config filename:"))
        filename_layout.addWidget(self.filename)
        self.outer_layout.addLayout(filename_layout)

        button_layout = QHBoxLayout()
        self.open_btn = QPushButton("Open Config")
        self.open_btn.clicked.connect(self.open_config)
        button_layout.addWidget(self.open_btn)
        self.save_btn = QPushButton("Save Config")
        self.save_btn.clicked.connect(self.save_config)
        button_layout.addWidget(self.save_btn)
        self.run_btn = QPushButton("Run Solver")
        self.run_btn.clicked.connect(self.run_solver)
        button_layout.addWidget(self.run_btn)
        self.outer_layout.addLayout(button_layout)

        # Connect radio buttons to enable/disable fields
        self.radio_calendar.toggled.connect(self.update_epoch_fields)
        self.update_epoch_fields()


        # for field in [self.moon_pa_file, self.ephemeris_file, self.gravfile, self.patch_point_file, self.initial_guess_from_file, self.object_name, self.leapseconds_file, self.mkspk_path, self.segment_id]:
        #     field.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

    # def plot_example(self):
    #     self.figure.clear()
    #     self.ax = self.figure.add_subplot(111)
    #     self.ax.clear()
    #     self.generate_defect_plot(self.ax, './defects_20251229220000_L2_S_NREVS=20.csv')
    #     self.ax.set_title("Defects")
    #     self.canvas.draw()

    def plot_example(self):

        try:
            self.plotter.remove_background_image()
        except Exception:
            pass  # Ignore if no background image exists

        self.plotter.clear()           # Remove all actors, charts, etc.
        # self.plotter.shape = (1, 1)    # Reset to a single subplot (or whatever shape you want)
        self.plotter.subplot(0, 0)     # Select the first (and only) subplot

        sys.path.insert(0, './python')
        from plot_utilities import get_pyvista_plotter
        get_pyvista_plotter('./traj_20251229220000_L2_S_NREVS=20.json', self.plotter, subplot_test=False)
        # # Example: plot a sphere
        # import pyvista as pv
        # sphere = pv.Sphere()
        # self.plotter.add_mesh(sphere, color='lightblue')
        self.plotter.reset_camera()
        self.plotter.update()   # <-- Ensure the widget updates
        self.plotter.render()   # <-- Force a redraw

    ###################################################################
    def generate_defect_plot(self, ax, filename : str, save : bool = True,
                             show : bool = False, ylabel : str = None,
                             showtitle : bool = True, dpi : int = 200):

        """
        Generate a plot of the constraint violations for a run.
        """

        idx = []
        r = []
        v = []
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            for i,row in enumerate(reader):
                if i==0: continue # skip header
                idx.append(float(i))
                r.append(float(row[0]))
                v.append(float(row[1]))

        ax.grid()
        ax.set_axisbelow(True)
        ax.set_xlabel("Segment interface number")
        if not ylabel:
            ylabel = 'Constraint violation magnitude'
        ax.set_ylabel(ylabel)
        if showtitle:
            ax.set_title(filename.strip('./'))
        ax.plot(idx,r,"b-",linewidth=2,label="r (km)")
        ax.plot(idx,v,"r-",linewidth=2,label="v (km/s)")
        ax.legend()

    def update_epoch_fields(self):
        calendar_enabled = self.radio_calendar.isChecked()
        for widget in [self.year, self.month, self.day, self.hour, self.minute, self.sec]:
            widget.setEnabled(calendar_enabled)
        self.et_ref.setEnabled(not calendar_enabled)

    def get_config(self):
        config = {
            "N_or_S": self.N_or_S.text(),
            "L1_or_L2": self.L1_or_L2.text(),
            "period": self.period.value(),
            "sa_lb": [self.sa_lb_0.value(), self.sa_lb_1.value()],
            "sa_ub": [self.sa_ub_0.value(), self.sa_ub_1.value()],
            "vm": [self.vm_0.value(), self.vm_1.value()],
            "sa_ns": self.sa_ns.value(),
            "sa_nt": self.sa_nt.value(),
            "sa_neps": self.sa_neps.value(),
            "n_revs": self.n_revs.value(),
            "fix_initial_time": self.fix_initial_time.isChecked(),
            "fix_initial_r": self.fix_initial_r.isChecked(),
            "fix_ry_at_end_of_rev": self.fix_ry_at_end_of_rev.value(),
            "fix_final_ry_and_vx": self.fix_final_ry_and_vx.isChecked(),
            "constrain_initial_rdot": self.constrain_initial_rdot.isChecked(),
            "fscale_rdot": self.fscale_rdot.value(),
            # "solver_mode": self.solver_mode.value(),
            "solver_mode": self.solver_mode_menu.currentData(),
            "rtol": self.rtol.text(),
            "atol": self.atol.text(),
            "nlesolver_tol": self.nlesolver_tol.text(),
            "use_splined_ephemeris": self.use_splined_ephemeris.isChecked(),
            "dt_spline_sec": self.dt_spline_sec.value(),
            "ephemeris_file": self.ephemeris_file.text(),
            "gravfile": self.gravfile.text(),
            "grav_n": self.grav_n.value(),
            "grav_m": self.grav_m.value(),
            # "grav_frame": self.grav_frame.value(),
            "grav_frame": self.grav_frame.currentData(),
            "moon_pa_file": self.moon_pa_file.text(),
            "patch_point_file": self.patch_point_file.text(),
            "patch_point_file_is_periapsis": self.patch_point_file_is_periapsis.isChecked(),
            "pointmass_central_body": self.pointmass_central_body.isChecked(),
            "include_pointmass_earth": self.include_pointmass_earth.isChecked(),
            "include_pointmass_sun": self.include_pointmass_sun.isChecked(),
            "include_pointmass_jupiter": self.include_pointmass_jupiter.isChecked(),
            "use_battin_gravity": self.use_battin_gravity.isChecked(),
            "initial_guess_from_file": self.initial_guess_from_file.text(),
            "generate_plots": self.generate_plots.isChecked(),
            "generate_trajectory_files": self.generate_trajectory_files.isChecked(),
            "generate_kernel": self.generate_kernel.isChecked(),
            "object_id": self.object_id.value(),
            "object_name": self.object_name.text(),
            "leapseconds_file": self.leapseconds_file.text(),
            "mkspk_path": self.mkspk_path.text(),
            "polynom_degree": self.polynom_degree.value(),
            "output_spk_type": self.output_spk_type.value(),
            "segment_id": self.segment_id.text(),
            "generate_guess_and_solution_files": self.generate_guess_and_solution_files.isChecked(),
            "generate_defect_file": self.generate_defect_file.isChecked(),
            "generate_json_trajectory_file": self.generate_json_trajectory_file.isChecked(),
            "generate_eclipse_files": self.generate_eclipse_files.isChecked(),
            "r_eclipse_bubble": self.r_eclipse_bubble.value(),
            "eclipse_dt_step": self.eclipse_dt_step.value(),
            "eclipse_filetype": self.eclipse_filetype.value(),
            "run_pyvista_script": self.run_pyvista_script.isChecked(),
            "generate_rp_ra_file": self.generate_rp_ra_file.isChecked()
        }
        if self.radio_calendar.isChecked():
            config.update({
                "epoch_mode": "calendar",
                "year": self.year.value(),
                "month": self.month.value(),
                "day": self.day.value(),
                "hour": self.hour.value(),
                "minute": self.minute.value(),
                "sec": self.sec.value()
            })
        else:
            config.update({
                "epoch_mode": "et",
                "et_ref": self.et_ref.value()
            })
        return config

    def set_config(self, config):
        # Mission tab
        self.N_or_S.setText(str(config.get("N_or_S", "N")))
        self.L1_or_L2.setText(str(config.get("L1_or_L2", "L1")))
        self.period.setValue(float(config.get("period", 0.0)))
        sa_lb = config.get("sa_lb", [0.0, 0.0])
        self.sa_lb_0.setValue(float(sa_lb[0]) if len(sa_lb) > 0 else 0.0)
        self.sa_lb_1.setValue(float(sa_lb[1]) if len(sa_lb) > 1 else 0.0)
        sa_ub = config.get("sa_ub", [0.0, 0.0])
        self.sa_ub_0.setValue(float(sa_ub[0]) if len(sa_ub) > 0 else 0.0)
        self.sa_ub_1.setValue(float(sa_ub[1]) if len(sa_ub) > 1 else 0.0)
        vm = config.get("vm", [0.0, 0.0])
        self.vm_0.setValue(float(vm[0]) if len(vm) > 0 else 0.0)
        self.vm_1.setValue(float(vm[1]) if len(vm) > 1 else 0.0)
        self.sa_ns.setValue(int(config.get("sa_ns", 0)))
        self.sa_nt.setValue(int(config.get("sa_nt", 0)))
        self.sa_neps.setValue(int(config.get("sa_neps", 0)))
        self.n_revs.setValue(int(config.get("n_revs", 0)))
        # Solver tab
        self.fix_initial_time.setChecked(bool(config.get("fix_initial_time", False)))
        self.fix_initial_r.setChecked(bool(config.get("fix_initial_r", False)))
        self.fix_ry_at_end_of_rev.setValue(int(config.get("fix_ry_at_end_of_rev", 0)))
        self.fix_final_ry_and_vx.setChecked(bool(config.get("fix_final_ry_and_vx", False)))
        self.constrain_initial_rdot.setChecked(bool(config.get("constrain_initial_rdot", False)))
        self.fscale_rdot.setValue(float(config.get("fscale_rdot", 0.0)))
        # self.solver_mode.setValue(int(config.get("solver_mode", 1)))
        solver_mode = int(config.get("solver_mode", 1))
        idx = next((i for i, (_, val) in enumerate(self.solver_mode_options) if val == solver_mode), 0)
        self.solver_mode_menu.setCurrentIndex(idx)

        # Integrator tab
        self.rtol.setText(str(config.get("rtol", "1.0e-12")))
        self.atol.setText(str(config.get("atol", "1.0e-12")))
        self.nlesolver_tol.setText(str(config.get("nlesolver_tol", "1.0e-6")))
        # Files tab
        self.use_splined_ephemeris.setChecked(bool(config.get("use_splined_ephemeris", False)))
        self.dt_spline_sec.setValue(int(config.get("dt_spline_sec", 3600)))
        self.ephemeris_file.setText(str(config.get("ephemeris_file", "")))
        self.gravfile.setText(str(config.get("gravfile", "")))
        self.grav_n.setValue(int(config.get("grav_n", 0)))
        self.grav_m.setValue(int(config.get("grav_m", 0)))
        # self.grav_frame.setValue(int(config.get("grav_frame", 0)))
        grav_frame_val = int(config.get("grav_frame", 0))
        idx = 0 if grav_frame_val == 1 else 1
        self.grav_frame.setCurrentIndex(idx)
        self.moon_pa_file.setText(str(config.get("moon_pa_file", "")))
        self.patch_point_file.setText(str(config.get("patch_point_file", "")))
        self.patch_point_file_is_periapsis.setChecked(bool(config.get("patch_point_file_is_periapsis", False)))
        self.pointmass_central_body.setChecked(bool(config.get("pointmass_central_body", False)))
        self.include_pointmass_earth.setChecked(bool(config.get("include_pointmass_earth", False)))
        self.include_pointmass_sun.setChecked(bool(config.get("include_pointmass_sun", False)))
        self.include_pointmass_jupiter.setChecked(bool(config.get("include_pointmass_jupiter", False)))
        self.use_battin_gravity.setChecked(bool(config.get("use_battin_gravity", False)))
        self.initial_guess_from_file.setText(str(config.get("initial_guess_from_file", "")))
        # Output tab
        self.generate_plots.setChecked(bool(config.get("generate_plots", False)))
        self.generate_trajectory_files.setChecked(bool(config.get("generate_trajectory_files", False)))
        self.generate_kernel.setChecked(bool(config.get("generate_kernel", False)))
        self.object_id.setValue(int(config.get("object_id", -50000)))
        self.object_name.setText(str(config.get("object_name", "")))
        self.leapseconds_file.setText(str(config.get("leapseconds_file", "")))
        self.mkspk_path.setText(str(config.get("mkspk_path", "")))
        self.polynom_degree.setValue(int(config.get("polynom_degree", 9)))
        self.output_spk_type.setValue(int(config.get("output_spk_type", 9)))
        self.segment_id.setText(str(config.get("segment_id", "")))
        self.generate_guess_and_solution_files.setChecked(bool(config.get("generate_guess_and_solution_files", False)))
        self.generate_defect_file.setChecked(bool(config.get("generate_defect_file", False)))
        self.generate_json_trajectory_file.setChecked(bool(config.get("generate_json_trajectory_file", False)))
        self.generate_eclipse_files.setChecked(bool(config.get("generate_eclipse_files", False)))
        self.r_eclipse_bubble.setValue(float(config.get("r_eclipse_bubble", 1.0)))
        self.eclipse_dt_step.setValue(float(config.get("eclipse_dt_step", 3600.0)))
        self.eclipse_filetype.setValue(int(config.get("eclipse_filetype", 2)))
        self.run_pyvista_script.setChecked(bool(config.get("run_pyvista_script", False)))
        self.generate_rp_ra_file.setChecked(bool(config.get("generate_rp_ra_file", False)))
        # Epoch tab
        epoch_mode = config.get("epoch_mode", "calendar")
        if epoch_mode == "calendar":
            self.radio_calendar.setChecked(True)
            self.year.setValue(int(config.get("year", 2025)))
            self.month.setValue(int(config.get("month", 1)))
            self.day.setValue(int(config.get("day", 1)))
            self.hour.setValue(int(config.get("hour", 0)))
            self.minute.setValue(int(config.get("minute", 0)))
            self.sec.setValue(float(config.get("sec", 0.0)))
        else:
            self.radio_et.setChecked(True)
            self.et_ref.setValue(float(config.get("et_ref", 0.0)))
        self.update_epoch_fields()

    def open_config(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Open Config File", "", "JSON Files (*.json);;All Files (*)")
        if fname:
            try:
                with open(fname, "r") as f:
                    config = json5.load(f)   # use json5 to allow comments
                self.set_config(config)
                self.filename.setText(fname)
                QMessageBox.information(self, "Success", f"Config loaded from {fname}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def save_config(self):
        config = self.get_config()
        fname = self.filename.text()
        try:
            with open(fname, "w") as f:
                json.dump(config, f, indent=2)
            QMessageBox.information(self, "Success", f"Config saved to {fname}")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def run_solver(self):
        fname = self.filename.text()
        self.solver_thread = SolverThread(["./run.sh", fname])
        self.progress_dialog = ProgressDialog(self)
        self.progress_dialog.cancel_btn.clicked.connect(self.cancel_solver)
        #self.solver_thread.output.connect(self.progress_dialog.console.appendPlainText)
        self.solver_thread.output.connect(lambda text: self.progress_dialog.console.appendPlainText(text.rstrip('\r\n')))
        self.solver_thread.finished.connect(self.solver_finished)
        self.solver_thread.start()
        self.progress_dialog.exec()

    def cancel_solver(self):
        if hasattr(self, 'solver_thread'):
            self.solver_thread.terminate_solver()
        self.progress_dialog.label.setText("Cancelling...")
        self.progress_dialog.cancel_btn.setEnabled(False)

    def solver_finished(self, returncode, stdout, stderr):
        self.plot_example()
        if hasattr(self, 'solver_thread'):
            self.solver_thread.terminate_solver()
        self.progress_dialog.label.setText("Done")
        QMessageBox.information(self, "HALO", "Solver finished")
        # self.progress_dialog.close()
        # if returncode == 0:
            # QMessageBox.information(self, "Solver Output", stdout)
        # else:
            # QMessageBox.critical(self, "Solver Error", stderr)


if __name__ == "__main__":

    #>.. i can't get this to work
    if sys.platform == "darwin":
        from AppKit import NSApplication, NSImage
        app_icon_path = "./media/logo.icns"
        app = NSApplication.sharedApplication()
        img = NSImage.alloc().initByReferencingFile_(app_icon_path)
        app.setApplicationIconImage_(img)

    app = QApplication(sys.argv)
    app.setApplicationName("Halo Solver")
    gui = HaloConfigGUI()
    gui.setWindowIcon(QIcon("./media/logo.png"))  # doesn't work
    gui.show()
    sys.exit(app.exec())