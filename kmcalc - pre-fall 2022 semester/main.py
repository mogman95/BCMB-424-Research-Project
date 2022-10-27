import auto_formatter
import analysis
import datetime as dt
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets, uic
from datetime import date
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
from uncertainties import ufloat
import sys
import os, platform
import global_settings
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
# import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
class Main(QtWidgets.QMainWindow):
    def __init__(self):
        super(Main, self).__init__()
        uic.loadUi(os.getcwd() + '/ui_files/kmmain-window.ui', self)
        self.setStyleSheet('QMainWindow{background-color: white;}')
        self.logo_label.setPixmap(QPixmap(u"" + os.getcwd() + "/ui_files/kmcalc_logo.png"))  # Sets the logo
        self.showMaximized()
        self.analyze_btn.clicked.connect(self.analysis)
        self.autoformat_btn.clicked.connect(self.autoformatter)
    def analysis(self):
        x = 'run analysis.py'
    def autoformatter(self):
        x = 'run auto_formatter.py'

app = QtWidgets.QApplication(sys.argv)  # Create an instance of QtWidgets.QApplication
window = Main()                # Create an instance of our class
app.exec_()                             # Start the application