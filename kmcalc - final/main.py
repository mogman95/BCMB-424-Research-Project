import auto_formatter
import analysis
from PyQt5.QtGui import *
from PyQt5 import QtWidgets, uic
import sys
import os, platform
import global_settings as gs
from PyQt5 import QtWidgets, QtCore, uic


class KmCalc(QtWidgets.QMainWindow):
    """
    Home window connecting the auto-formatting and analysis programs.
    """
    def __init__(self):
        """
        Shows and formates the home window and initializes the subprogram windows.
        """
        super(KmCalc, self).__init__()
        uic.loadUi(gs.uidir + 'home-window.ui', self)
        self.setStyleSheet('QMainWindow{background-color: white;}')
        self.logo_label.setPixmap(QPixmap(gs.uidir + "kmcalc_logo2.png").scaled(750, 750, aspectRatioMode=QtCore.Qt.KeepAspectRatio,transformMode=QtCore.Qt.SmoothTransformation))  # Sets the logo
        self.mwfg = self.frameGeometry() # Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  # Center window
        self.mwfg.moveCenter(self.cp)  # Center window
        self.move(self.mwfg.topLeft())  # Center window


        gs.rows_errDialog = auto_formatter.rows_errDialog()
        gs.DataConf = auto_formatter.DataConf()
        gs.colMap = auto_formatter.colMap()
        gs.DataSelection = auto_formatter.DataSelection(os.getcwd())
        gs.main_window = analysis.MainWindow(gs.appdir)
        gs.results_window = analysis.ResultsWindow()

        self.autoformat_btn.clicked.connect(self.autoformatter) # Runs the auto-formatting program when its respective button is pressed
        self.analyze_btn.clicked.connect(self.analysis) # Runs the analysis program when its respective button is pressed


        self.show() # Show the window

    def analysis(self):
        """
        Shows the analysis main window and establishes return home button connections.
        """
        gs.main_window.show()
        gs.results_window.rtn_home_btn.clicked.connect(self.rtn_home)
        gs.main_window.return_home_btn.clicked.connect(self.rtn_home)

        self.hide() # Hide main

    def autoformatter(self):
        """
        Shows the auto-formatter data selection window and establishes return button connections.
        """
        gs.DataSelection.show()
        gs.colMap.rtn_home_btn.clicked.connect(self.rtn_home)
        gs.DataSelection.return_home_btn.clicked.connect(self.rtn_home)

        self.hide() # Hide main
    
    def rtn_home(self):
        """
        Hides the window on which the return home button was pressed and shows the home window.
        """
        if gs.colMap.isVisible():
            gs.colMap.hide()
            self.show()
        if gs.DataSelection.isVisible():
            gs.DataSelection.hide()
            self.show()
        if gs.main_window.isVisible():
            gs.main_window.hide()
            self.show()
        if gs.results_window.isVisible():
            gs.results_window.go_back(return_home=True)
            self.show()

        


def main():
    """
    Executes the main program.
    """
    if hasattr(sys, 'frozen'):
        gs.appdir = sys.executable
        gs.uidir = sys.executable

        if platform.system() == 'Windows':
            gs.appdir = gs.appdir[:gs.appdir.rfind("\\") + 1]
            gs.uidir += "\\ui_files\\"
        else:
            gs.appdir = gs.appdir[:gs.appdir.rfind/("/") + 1]
            gs.uidir += sys.executable + "/ui_files/"
    else:
        gs.appdir = os.path.dirname(os.path.abspath(__file__))
        gs.uidir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            gs.appdir += '\\'
            gs.uidir += "\\ui_files\\"
        else:
            gs.appdir += '/'
            gs.uidir += "/ui_files/"

    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = QtWidgets.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("KmCalc")
    print("Current app directory is: " + gs.appdir + "\n")
    print("Current ui directory is: " + gs.uidir + "\n")
    window = KmCalc()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
