from datetime import date
# from matplotlib.pyplot import cm, figure
import pandas as pd
import numpy as np
from scipy.stats import linregress
# from statistics import mean, stdev
from scipy.optimize import curve_fit
from uncertainties import ufloat
import sys
import os, platform
import global_settings
# import io
from PyQt5 import QtWidgets, Qt, QtGui, QtCore, uic
# import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar

# =========================================================================================
# CLASS NAME: MainWindow
# Inputs: Takes in the path information of the program
# Outputs: KmCalc job (output to files and graphs)
# =========================================================================================


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self,appdir):
        super(MainWindow, self).__init__()
        uic.loadUi(global_settings.appdir + "/ui_files/kmcalc_main.ui", self)
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.mwfg.moveCenter(self.cp)  ##Center window
        self.move(self.mwfg.topLeft())  ##Center window
        self.fill_boxes()

        # --- Variable Initialization --- #
        self.input_path = ""
        self.output_path = ""
        self.linear_tol = float()
        self.verbose_level = True

        # --- Style Modifications --- #
        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(94,198,201);
                        border-radius: 5px;
                        font: 14pt "Arial";
                        font: bold;
                        margin-top: 10px;}"""

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2").replace("rgb(94,198,201)","rgb(240,188,193)"))

        # --- Button Modifications --- /Applications/Qt Designer.app/Contents/MacOS/Designer#
        self.close_button.clicked.connect(self.close)
        self.browse_button.clicked.connect(self.browse_file)
        self.submit_button.clicked.connect(self.submit)

        # --- Show GUI (should be last command in init) --- #
        self.show()

    def gather_args(self):
        if self.output_edit.text() == "":
            self.output_path = self.input_path.strip(".xlsx")[0] + "-output.xlsx"
        elif (self.output_edit.text().endswith(".xlsx") or self.output_edit.text().endswith(".csv")):
            print('line 67',global_settings.appdir)
            self.output_path = global_settings.appdir + self.output_edit.text()
        else:
            self.output_path = global_settings.appdir + self.output_edit.text() + ".xlsx"
        if self.verbose_box.currentText() == "Yes":
            self.verbose_level = True
        else:
            self.verbose_level = False

        self.linear_tol = float(self.tolerance_box.value())
        if self.linear_tol == "": ### <= user should have to choose for transparency and consistency (and/or have option to autofind lin. tol.)
            self.linear_tol = float(1.5)

    def submit(self):
        if self.input_path != '':
            self.gather_args()
            self.find_parameters(self.input_path, self.output_path, self.verbose_level, self.linear_tol)
        else:
            self.input_edit.setText('Please select a file to analyze')
    
    def fill_boxes(self):
        yes_no_list = ["Yes","No"]
        for item in yes_no_list:
            self.plot_box.addItem(item)
            self.verbose_box.addItem(item)

    def browse_file(self):
        # self.input_path = "/home/ddooley/software/KmCalc/examples/1/ex1.xlsx"
#         self.input_path = "/Users/ddooley/bioinformatics_packages/individual_packages/kmcalc/examples/3/ex3.xlsx"
        filed = QtWidgets.QFileDialog()
        myFile = QtWidgets.QFileDialog.getOpenFileName(filed, caption="Please select an input file...", # Sets the title of the window
                                                directory=os.getcwd(),        # Sets the starting directory as the current working directory
                                                filter='excel (*.xlsx)')           # Sets the file filter
        if (myFile[0] != ""):
            my_file = str(myFile[0])
            self.input_edit.setText(my_file)
            self.input_path = my_file
            global_settings.input_path = my_file
            if platform.system() == "Windows":
                global_settings.input_path = global_settings.input_path.replace("/","\\")
            else:
                global_settings.input_path = global_settings.input_path.replace("\\","/")

    """ THIS SECTION HOLDS THE ACTUAL FUNCTIONS THAT CALCULATE KINETIC PARAMETERS """

    # Main function
    def find_parameters(self,input_path, output_path, verbose_level, linear_region_tol):
        global verbose_bool 
        self.hide() ### Hide main window
        global_settings.loading_window.show() ### Show loading window
        global_settings.loading_window.progressBar.setValue(0) ### Update progress Bar
        
        QtCore.QCoreApplication.processEvents()
        self.print_citation() ### Print citation information to output window

        verbose_bool = verbose_level
        raw_df, col_map, conversions = self.read_input(input_path)
        global_settings.loading_window.progressBar.setValue(10) ### Update progress Bar
        QtCore.QCoreApplication.processEvents()

        ### Find units
        for i, col in enumerate(raw_df.columns):
            item = str(col).replace(' ', '')
            if 'time' in str(item).lower():
                self.time_units = item[5:len(item) - 1]  # Assumes 'Time (units)' format
                break
        ### The conc units should be the same but I included them
        # so that potential error checking could be implimented in the future.
        # One variable (e.g. conc_units) could be used instead
        self.sub_conc_units = col_map.columns[3][13:15]
        self.enz_conc_units = col_map.columns[4][10:12]
        print(self.time_units, self.sub_conc_units, self.enz_conc_units)
        global_settings.loading_window.progressBar.setValue(30)  ### Update progress Bar
        QtCore.QCoreApplication.processEvents()

        df = raw_df[[f'Time ({self.time_units})']+list(col_map['Well'])] # Only keep columns noted in column map #==> Needs way to find time units
        rates_abs_min, rates_r2, linear_region = self.find_rates(df, linear_region_tol, col_map)
        global_settings.loading_window.progressBar.setValue(40) ### Update progress Bar
        QtCore.QCoreApplication.processEvents()

        # Convert rates to target units
        # print(rates_abs_min)
        rates_min = rates_abs_min*conversions.loc['absorbance2concentration', 'value']
        # print(rates_min)
        rates = rates_min*conversions.loc['time_conversion', 'value']
        # print(rates)
        global_settings.loading_window.progressBar.setValue(60) ### Update progress Bar
        QtCore.QCoreApplication.processEvents()

        # format replicates as: # {sample_id: {rep1: conc_rate_df, rep2:...}}
        replicates = self.split_replicates_and_subtract_blank(rates, col_map)
        out_df, out_dict = self.get_output_df(replicates, conversions)
        global_settings.loading_window.output_window.append(output_path)
        rates_out = self.write_report(output_path, out_df, replicates, rates_r2, col_map, linear_region)
        global_settings.loading_window.progressBar.setValue(80) ### Update progress Bar
        QtCore.QCoreApplication.processEvents()

        global_settings.results_window.generate_plots(output_path, rates_out, out_dict)
        global_settings.loading_window.progressBar.setValue(90) ### Update progress Bar
        QtCore.QCoreApplication.processEvents()
        
        global_settings.loading_window.progressBar.setValue(100) ### Update progress Bar
        QtCore.QCoreApplication.processEvents()
        global_settings.loading_window.continue_fun()

    def read_input(self,input_path):
        if input_path.endswith('xlsx'):
            raw_df = pd.read_excel(input_path, sheet_name='Raw Data')
            col_map = pd.read_excel(input_path, sheet_name='Column Map')
            conversions = pd.read_excel(input_path, sheet_name='Conversion Factors', index_col='name')
        else:
            raw_df = pd.read_csv(os.path.join(input_path, 'raw_data.csv'))
            col_map = pd.read_csv(os.path.join(input_path, 'column_map.csv'))
            conversions = pd.read_csv(os.path.join(input_path, 'conversion_factors.csv'), index_col='name')
        return raw_df, col_map, conversions

    #def get_output_path(output_path, input_path):
    #    if output_path in [None, 'default']:
    #        filename, file_extension = os.path.splitext(input_path)
    #        if input_path.endswith('xlsx'):
    #            output_path = '{}_output.xlsx'.format(filename)
    #        else:
    #            output_path = os.path.dirname(input_path[:-1])
    #    return output_path

    def write_report(self,output_path, out_df, replicates, rates_r2, col_map, linear_region):

        #Use the rates from replicates without blank:
        frames = []
        for k,v in replicates.items():
            for k1,v1 in v.items():
                frames.append(v1)
        rates_r = pd.concat(frames)['rate'].round(4)
        rates_r2_r = rates_r2.round(4)
        rates_r.name =f'rate({self.sub_conc_units}/s)'
        rates_out = col_map.set_index('Well').join(rates_r).join(rates_r2_r).join(linear_region)

        if output_path.endswith('xlsx'):
            with pd.ExcelWriter(output_path, engine='xlsxwriter') as w:
                out_df.to_excel(w, sheet_name='parameters')
                rates_out.to_excel(w, sheet_name='rates')
                w.save()
        else:
                out_df.to_csv(os.path.join(output_path,'parameters.csv'))
                rates_out.to_csv(os.path.join(output_path,'rates.csv'))
        if verbose_bool:
            print('Output written to: {}'.format(os.path.abspath(output_path)))


        return rates_out

    
    def print_citation(self):
        global_settings.loading_window.output_window.append('If you use any of this software, please cite:')
        global_settings.loading_window.output_window.append('CITATION_PLACEHOLDER')
        QtCore.QCoreApplication.processEvents()

    def find_rates(self,df, linear_region_tol, col_map):
        # todo: make this method return a frame instead of series
        time = df[f'Time ({self.time_units})']
        rates = pd.Series(index=df.columns[1:], name='rate')
        rates_r2 = pd.Series(index=df.columns[1:], name='R2')
        linear_region = pd.Series(index=df.columns[1:], name=f'linear_region({self.time_units})')
        if ("start_linear_region" in col_map.columns) and ("end_linear_region" in col_map.columns): # Linear region is specified by the user
            def get_linear_region(time, col, linear_region_tol, col_map):
                return self.read_linear_region(time,col, col_map)
        else:
            def get_linear_region(time, col, linear_region_tol, col_map):
                return self.find_linear_region(time, col, linear_region_tol)

        for col_id in df:
            col = df[col_id]
            if col.name == f'Time ({self.time_units})':
                continue
            in_linear_region = get_linear_region(time, col, linear_region_tol, col_map)
            result = linregress(x=time[in_linear_region], y=col[in_linear_region])
            rates[col_id] = result[0]
            rates_r2[col_id] = result[2]
            linear_region[col_id] = '{}-{}'.format(time[in_linear_region].iloc[0],time[in_linear_region].iloc[-1])
        return rates, rates_r2, linear_region

    def read_linear_region(self,time, col, col_map):
        start = col_map.loc[ col_map['Well'] == col.name, 'start_linear_region'].values[0]
        end = col_map.loc[ col_map['Well'] == col.name, 'end_linear_region'].values[0]
        return (time > start) & (time < end)

    def find_linear_region(self,time, col, linear_region_tol):# create theoretical data
        # Convolution filter adapted from stackexchange
        # Alternatively one can fit a spline: https://stackoverflow.com/questions/40226357/second-derivative-in-python-scipy-numpy-pandas
        #

        y_n = col
        x = time

        # create convolution kernel for calculating  the smoothed second order derivative
        smooth_width = 59
        x1 = np.linspace(-3,3,smooth_width)
        #norm = np.sum(np.exp(-x1**2)) * (x1[1]-x1[0]) # ad hoc normalization
        y1 = (4*x1**2 - 2) * np.exp(-x1**2) / smooth_width *8#norm*(x1[1]-x1[0])

        # calculate second order deriv.
        y_conv = np.convolve(y_n, y1, mode="same")

        # Identify linear regions
        is_linear = abs(y_conv) <= linear_region_tol

        # Note: Algorithm below may need adjustment if the first linear region detected is very small.
        # Select first contiguous linear region:
        in_linear_region = np.full((len(is_linear),1), False, dtype=bool)
        first = False
        for idx, item in enumerate(is_linear):
            if item:
                first = True
            else:
                if first: # item is not linear but first linear region was already mapped
                    break
            if first: # we are in the linear region
                in_linear_region[idx] = True
        in_linear_region = pd.DataFrame(in_linear_region)[0] # make a series for boolean indexing
        return in_linear_region

    def split_replicates_and_subtract_blank(self,rates, col_map):
        # format replicates as: # {sample_id: {rep1: conc_rate_df, rep2:...}}
        replicates = {}
        for sample_id in col_map['Sample ID'].unique():
            enz_map = col_map[col_map['Sample ID'] == sample_id]
            rep_dict = {}
            for rep in enz_map['Replicate'].unique():
                rep_map = enz_map[enz_map['Replicate'] == rep]
                rep_df = rep_map.set_index('Well').join(rates)
                # subtract blank
                blank_value = rep_df[rep_df[f'[Substrate] ({self.sub_conc_units})'] == 0]['rate'].values
                rep_df['rate'] = rep_df['rate'] - blank_value
                rep_dict[rep] = rep_df
            replicates[sample_id] = rep_dict
        return replicates


    def get_output_df(self,replicates, conversions):
        n_output_decimals = 4
        out_dict = {}
        for sample_id, reps in replicates.items():
            # From a statistical point of view it makes more sense to fit the model to all available data instead of
            # treating replicates individually, thus:
            all_rep_df = pd.concat([rep_df for id, rep_df in reps.items()])
            km, vmax = self.find_km_vmax(all_rep_df)
            # Enzyme conc (ug/mL
            kcat = vmax/all_rep_df.loc[
                all_rep_df['Sample ID'].str.contains(sample_id).index[0], # .astype(str)
                f'[Enzyme] ({self.enz_conc_units})']
            kcatkm = kcat/km
            out_dict[sample_id] = {f'km ({self.sub_conc_units})':km, f'vmax ({self.sub_conc_units}/s)':vmax, f'kcat (1/s)':kcat, f'kcat/km (1/s/{self.sub_conc_units})':kcatkm}
        out_df = pd.DataFrame(out_dict).transpose()
        fstr = '{:.' + str(n_output_decimals) + 'f}'
        for i in range(out_df.shape[0]): # looping migh not be ideal by df.apply did not work well
            for j in range(out_df.shape[1]):
                out_df.iloc[i,j] = fstr.format(out_df.iloc[i,j])

        return out_df, out_dict


    def mm_model(self,s, vmax, km):
        return (vmax*s)/(km + s)


    def find_km_vmax(self,rep_df):
        # Use robust non-linear regression
        # The paremeter covariance estimate procedure is explained here: https://stackoverflow.com/questions/14854339/in-scipy-how-and-why-does-curve-fit-calculate-the-covariance-of-the-parameter-es
        # bounds and inital guesses on km and vmax parameters (concentration units), from Cell Biology by the Numbers: http://book.bionumbers.org/how-many-reactions-do-enzymes-carry-out-each-second/
        bounds = (np.full([2,], 0), np.full([2,], 10e4)) # WARNING: THIS BOUNDS ASSUME CERTAIN INPUT UNITS!
        x0 = np.full([2,], 0.13) # vmax is not given in the refernce but rather kcat

        popt, pcov = curve_fit(self.mm_model, xdata=rep_df[f'[Substrate] ({self.sub_conc_units})'],
                            ydata=rep_df['rate'], p0=x0, bounds=bounds, loss='soft_l1')
        perr = np.sqrt(np.diag(pcov)) # Standard deviation of the parameters.

        vmax = ufloat(popt[0], perr[0])
        km = ufloat(popt[1], perr[1])

        return km, vmax


"""
    # this function prepares everything for the generate library function
    # it is very similar to the gather settings, how ever it stores the data instead of calling the Annotation Window class
    # it moves the data onto the generateLib function, and then opens that window
    def fill_annotation_dropdown(self):
        temp_list = list()
        self.annotation_files.clear()
        for file in os.listdir(global_settings.CSPR_DB):
            if ".gbff" in file:
                temp_list.append(str(file))
        temp_list.sort(key=str.lower)
        for file in temp_list:
            self.annotation_files.addItem(file)
"""

#    def open_ncbi_web_page(self):
#        webbrowser.open('https://www.ncbi.nlm.nih.gov/', new=2)
#
#
#    def open_casper2_web_page(self):
#        webbrowser.open('http://casper2.org/', new=2)
#
#
#    def visit_repo_func(self):
#        webbrowser.open('https://github.com/TrinhLab/CASPERapp')

    # this if the function that is called when the user closes the program entirely.
    # so far I only know of 4 spots that can do this
    #       1. mainWindow
    #       2. annotationsWindow
    #       3. Results
    #       4. Multitargetting

class LoadingWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(LoadingWindow, self).__init__()
        uic.loadUi(global_settings.appdir + "/ui_files/kmloading_window.ui", self)
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.mwfg.moveCenter(self.cp)  ##Center window
        self.move(self.mwfg.topLeft())  ##Center window
        self.back_button.clicked.connect(self.go_back)
        self.continue_button.clicked.connect(self.continue_fun)

        # --- Variable Initialization --- #

        # --- Style Modifications --- #
        # groupbox_style = """
        # QGroupBox:title{subcontrol-origin: margin;
        #                 left: 10px;
        #                 padding: 0 5px 0 5px;}
        # QGroupBox#Step1{border: 2px solid rgb(94,198,201);
        #                 border-radius: 5px;
        #                 font: 14pt "Arial";
        #                 font: bold;
        #                 margin-top: 10px;}"""

        # self.Step1.setStyleSheet(groupbox_style)
        # self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2").replace("rgb(94,198,201)","rgb(240,188,193)"))

        # --- Button Modifications --- #
        # self.back_button.clicked.connect(INSERT)
        # self.continue_button.clicked.connect(INSERT)

        # --- Show GUI (should be last command in init) --- #
        self.hide()
    
    def go_back(self):
        self.hide()
        self.output_window.clear()
        global_settings.main_window.show()

    def continue_fun(self):
        self.hide()
        self.output_window.clear()
        #global_settings.results_window.populate_graphs()
        QtCore.QCoreApplication.processEvents()
        global_settings.results_window.show()



class ResultsWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(ResultsWindow, self).__init__()
        uic.loadUi(global_settings.appdir + "/ui_files/kmresults_window.ui", self)
        self.name_window = name_window()
        self.plots_generated = False
        self.mwfg = self.frameGeometry()  ##Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  ##Center window
        self.mwfg.moveCenter(self.cp)  ##Center window
        self.move(self.mwfg.topLeft())  ##Center window

        ### Initialize layouts for figures
        self.widget_layout = QtWidgets.QVBoxLayout()
        self.widget_layout.setContentsMargins(0,0,0,0)
        self.widget_2_layout = QtWidgets.QVBoxLayout()
        self.widget_2_layout.setContentsMargins(0,0,0,0)

        ### Connect functions to buttons
        self.back_button.clicked.connect(self.go_back)
        self.save_all.triggered.connect(self.save_all_figures)

        groupbox_style = """
                QGroupBox:title{subcontrol-origin: margin;
                                left: 10px;
                                padding: 0 5px 0 5px;}
                QGroupBox#all_box{border: 2px solid rgb(94,198,201);
                                border-radius: 5px;
                                font: 14pt "Arial";
                                font: bold;
                                margin-top: 10px;}"""

        self.all_box.setStyleSheet(groupbox_style)
        self.individual_box.setStyleSheet(groupbox_style.replace("all_box", "individual_box").replace("rgb(94,198,201)","rgb(240,188,193)"))


    def go_back(self):
        self.hide()
        self.plots_generated = False
        self.tabWidget.clear()
        global_settings.main_window.show()

    def generate_plots(self,output_path, rates_out, out_dict):
        self.sample_ids = []
        self.enz_ids = []
        self.sub_ids = []
        self.canvas_list = []

        ### Clear out old widgets in layout
        self.widget_canvas = MplCanvas(self, width=5, height=4, dpi=100) # Initialize new Canvas
        toolbar = NavigationToolbar(self.widget_canvas, self) # Initialize toolbar for widget 1
        self.widget_layout.addWidget(toolbar) # Add canvas to widget 1 layout
        self.widget_layout.addWidget(self.widget_canvas) # Add canvas to widget 1 layout
        self.widget.setLayout(self.widget_layout) # Add canvas (with layout) to widget 1 object

        if output_path.endswith('xlsx'):
            output_dir = os.path.join(os.path.dirname(output_path), 'plots')
        else:
            output_dir = os.path.join(output_path, 'plots')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        counter = 0
        plot_styles = ["o","v","^","s","D","P","8","H","X","*","p"]
        colors = ["cornflowerblue","indianred","gold","seagreen","rebeccapurple","orangered","lightslategrey","cadetblue","steelblue","navy","dimgray"]
        for sample_id, params in out_dict.items():
            self.sample_ids.append(sample_id) # Add enzyme ID to list
            sepi = sample_id.find('~') # Separator index: finds the location of the ~ separator
            self.enz_ids.append(sample_id[0:sepi])
            self.sub_ids.append(sample_id[sepi+1:len(sample_id)])
            # Initialize widget, canvas, layout, and toolbar for each new sample
            temp_widget = QtWidgets.QWidget() # Create new widget
            temp_canvas = MplCanvas(self, width=5, height=4, dpi=100) # Initialize new Canvas
            temp_layout = QtWidgets.QVBoxLayout() # Create new layout for widget
            temp_layout.setContentsMargins(0,0,0,0)
            temp_toolbar = NavigationToolbar(temp_canvas, self) # Initialize toolbar for new graph 
            temp_layout.addWidget(temp_toolbar) # Add toolbar to widget layout
            temp_layout.addWidget(temp_canvas) # Add canvas to widget layout
            temp_widget.setLayout(temp_layout)

            erates = rates_out[rates_out['Sample ID'] == sample_id]
            # Generate figures
            for i, rep in enumerate(erates['Replicate'].unique()):
                rerates = erates[erates['Replicate']==rep]
                # Get plot variables
                self.e_rates = rerates[f'rate({global_settings.main_window.sub_conc_units}/s)']
                self.s_conc = rerates[f'[Substrate] ({global_settings.main_window.sub_conc_units})']
                temp_canvas.axes.scatter(self.s_conc, self.e_rates,c=colors[i],label='Replicate {}'.format(rep)) ### Add each replicate to individual and total plots
                self.widget_canvas.axes.scatter(self.s_conc, self.e_rates,marker=plot_styles[counter],c=colors[counter]) ### Add each replicate to individual and total plots

            self.ps_conc = np.linspace(min(self.s_conc), max(self.s_conc), 200)
            km = params[f'km ({global_settings.main_window.sub_conc_units})']
            vmax = params[f'vmax ({global_settings.main_window.sub_conc_units}/s)']
            self.pe_rates = [self.mm_model(x, km=km.nominal_value, vmax=vmax.nominal_value) for x in self.ps_conc]
            temp_canvas.axes.plot(self.ps_conc, self.pe_rates, c="k", ls='--',
                     label='$\\bf{{Model:}}$\nK$_m$={:.3fP}\nV$_{{max}}$={:.3fP}'.format(km,vmax))

            # Plot enzyme's data on total figure
            self.widget_canvas.axes.plot(self.ps_conc, self.pe_rates, c=colors[counter], ls='--',
            label=(r'$\bf{}:$'.format("{" + sample_id.replace("_","-") + "}")+'\nK$_m$={:.3fP}\nV$_{{max}}$={:.3fP}'.format(km,vmax)))

#                     label=r'$\bf{{sample_id}}$' + ' Model:\n km={:.1f}\n vmax={:.4f}'.format(km,vmax))

            # Individual plot details
#            temp_canvas.axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
            temp_canvas.axes.legend(fontsize=8).set_draggable(state=True)
            temp_canvas.axes.set_xlabel(f'Initial substrate ({global_settings.main_window.sub_conc_units})')
            temp_canvas.axes.set_ylabel(f'Rate ({global_settings.main_window.sub_conc_units}/s)')

            # Total plot details
#            self.widget_canvas.axes.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=8)
            self.widget_canvas.axes.legend(fontsize=8).set_draggable(state=True)
            self.widget_canvas.axes.set_xlabel(f'Initial substrate ({global_settings.main_window.sub_conc_units})')
            self.widget_canvas.axes.set_ylabel(f'Rate ({global_settings.main_window.sub_conc_units}/s)')

            # find max rates across ALL enzymes to make all plots on the same scale
            self.widget_canvas.axes.set_ylim([0, max(rates_out[f'rate({global_settings.main_window.sub_conc_units}/s)'])*1.1])
            self.tabWidget.addTab(temp_widget,sample_id) # Add new tab and append tab_canvas to it

#            temp_canvas.axes.legend().set_draggable(state=True) # Make legends draggable
            self.canvas_list.append(temp_canvas)

            counter +=1 # Iterate counter

        self.plots_generated = True # Set bool for whether plots have been generated to true.

    def mm_model(self,s, vmax, km):
        return (vmax*s)/(km + s)

    def save_all_figures(self):
        if self.plots_generated:
            today = date.today()
            figure_count = self.tabWidget.count() ### Number of individual figures in tab widget 
            tmp = "all_enzymes_"+str(today.strftime("%m/%d/%y")).replace("/","_")+".png"
            self.name_window.rename_table.setRowCount(1+figure_count)

            ### Add table row for the graph with all enzymes
            item = QtWidgets.QTableWidgetItem(tmp)
            item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.name_window.rename_table.setItem(0, 0, item)
            self.name_window.rename_table.setCellWidget(0, 1, QtWidgets.QLineEdit())

            cnt = 1
            for i in range(figure_count):
                item = QtWidgets.QTableWidgetItem(self.tabWidget.tabText(i)+tmp+".png")
                item.setFlags(QtCore.Qt.ItemIsEnabled)
                self.name_window.rename_table.setItem(cnt, 0, item)
                self.name_window.rename_table.setCellWidget(cnt, 1, QtWidgets.QLineEdit())
                cnt += 1

            header = self.name_window.rename_table.horizontalHeader()
            self.name_window.rename_table.resizeColumnsToContents()
            header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
            header.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
            self.name_window.resize(self.name_window.sizeHint())

            ### Open dialog window to let them save figures
            msgBox = QtWidgets.QMessageBox()
            msgBox.setIcon(QtWidgets.QMessageBox.Icon.Question)
            msgBox.setWindowTitle("Name Figures:")
            msgBox.setText("Would you like to rename the figures?")
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.Yes)
            msgBox.addButton(QtWidgets.QMessageBox.StandardButton.No)
            msgBox.exec()

            if msgBox.result() == QtWidgets.QMessageBox.Yes:
                self.current_names = False
                self.name_window.show()
            else:
                print("Saving figures...")
                self.current_names = True
                self.name_window.save_figures()
        else:
            print("Plots have not been generated")
            

#     def populate_graphs(self):
#         ### Clear out old widgets in layout
#         self.widget_canvas = MplCanvas(self, width=5, height=4, dpi=100) # Initialize new Canvas
#         self.widget_2_canvas = MplCanvas(self, width=5, height=4, dpi=100) # Initialize new Canvas
#         toolbar = NavigationToolbar(self.widget_canvas, self) # Initialize toolbar for widget 1
#         toolbar_2 = NavigationToolbar(self.widget_2_canvas, self) # Initialize toolbar for widget 1
#         self.widget_layout.addWidget(toolbar) # Add canvas to widget 1 layout
#         self.widget_layout.addWidget(self.widget_canvas) # Add canvas to widget 1 layout
#         self.widget_2_layout.addWidget(toolbar_2) # Add canvas to widget 2 layout
#         self.widget_2_layout.addWidget(self.widget_2_canvas) # Add canvas to widget 2 layout
#         self.widget.setLayout(self.widget_layout) # Add canvas (with layout) to widget 1 object
#         self.widget_2.setLayout(self.widget_2_layout) # Add canvas (with layout) to widget 2 object

#         """ Populate widget 1 first """
# #        self.widget_canvas.axes.scatter(global_settings.main_window.ps_conc, global_settings.main_window.pe_rates, c='k', ls='--',
#         ### The following statements are plottings / formatting for the graph
#         rep = ["1","2","3"]

#         self.widget_canvas.axes.scatter(global_settings.main_window.s_conc, global_settings.main_window.e_rates, label='Replicate {}'.format(rep))
# #        self.widget_canvas.axes.set_xticklabels(x_labels)
#         self.widget_canvas.axes.set_xlabel('Test', fontsize=10)
#         self.widget_canvas.axes.set_ylabel('Test', fontsize=10)
#         self.widget_canvas.axes.set_title('Test',fontsize=10)
#         self.widget_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
#         self.widget_canvas.draw()

#         # self.seed_canvas.axes.yaxis.set_major_locator(MaxNLocator(integer=True))
#         # self.seed_canvas.axes.set_ylim(0, max(y) + 1)
#         # self.seed_canvas.axes.set_xticks(x)
#         # self.seed_canvas.axes.set_xticklabels(x_labels)
#         # if len(x_labels) > 10:
#         #     tick_spacing = round(len(x_labels)/10)
#         #     for i, t in enumerate(self.seed_canvas.axes.get_xticklabels()):
#         #         if (i % tick_spacing) != 0:
#         #             t.set_visible(False)
#         # self.seed_canvas.axes.set_xlabel('Chromosome', fontsize = 10)
#         # self.seed_canvas.axes.set_ylabel('Number of Repeats', fontsize=10)
#         # self.line_canvas.axes.set_title('Repeats per Scaffold/Chromosome',fontsize=10)
#         # self.line_canvas.axes.tick_params(axis='both', which='major', labelsize=8)
#         # self.line_canvas.draw()

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi,tight_layout=True)
        self.axes = fig.add_subplot(111)
        self.axes.clear()
        super(MplCanvas, self).__init__(fig)

class name_window(QtWidgets.QMainWindow):
    def __init__(self):
        super(name_window, self).__init__()
        uic.loadUi(global_settings.appdir + "/ui_files/kmname_window.ui", self)
        self.setWindowTitle("Name Figures")
        self.rename_table.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.rename_table.setColumnCount(2)
        self.rename_table.setHorizontalHeaderLabels(['Current Name', 'New Name'])

        ### Set pixel width for scroll bars
        self.rename_table.verticalScrollBar().setStyleSheet("width: 16px;")
        self.rename_table.horizontalScrollBar().setStyleSheet("height: 16px;")

        ### Connect functions to buttons
        self.submit_button.clicked.connect(self.save_figures)
        self.go_back_button.clicked.connect(self.go_back)

        self.label.setVisible(False) # Hide "Saving Figures:" label until Submit is clicked
        self.hide()
    
    def save_figures(self):
        pb = self.progressBar
        output_dir = os.path.join(global_settings.appdir, 'plots')
        self.label.setVisible(True)
        pb.setValue(20)
        increment = int(round(80/(len(global_settings.results_window.canvas_list))))
        QtCore.QCoreApplication.processEvents()
        
        if global_settings.results_window.current_names: # If user chose not to rename figures
            """ Start by saving the all_enzymes figure """
            plot_name = self.rename_table.item(0,0).text()
            global_settings.results_window.widget_canvas.axes.figure.savefig(os.path.join(output_dir,plot_name), bbox_inches = 'tight', dpi=500)

            """ Now save all of the individual enzyme figures """
            for i, canvas in enumerate(global_settings.results_window.canvas_list):
                i+=1
                plot_name = self.rename_table.item(i,0).text()
                canvas.axes.figure.savefig(os.path.join(output_dir,plot_name), bbox_inches = 'tight', dpi=500)

            print("Finished saving figures.")
        else: # If user chose to rename figures
            """ Start by saving the all_enzymes figure """
            if self.rename_table.cellWidget(0,1).text() != "": # If rename entry is not empty, save under the entered name
                plot_name = self.rename_table.cellWidget(0,1).text()
                global_settings.results_window.widget_canvas.axes.figure.savefig(os.path.join(output_dir,plot_name), bbox_inches = 'tight', dpi=500)
                pb.setValue(pb.value()+increment) # Update progress bar
                QtCore.QCoreApplication.processEvents()
            else: # If rename entry is empty, save under original name
                plot_name = self.rename_table.item(0,0)
                global_settings.results_window.widget_canvas.axes.figure.savefig(os.path.join(output_dir,plot_name), bbox_inches = 'tight', dpi=500)
                pb.setValue(pb.value()+increment) # Update progress bar
                QtCore.QCoreApplication.processEvents()

            """ Now save all of the individual enzyme figures """
            for i, canvas in enumerate(global_settings.results_window.canvas_list):
                i+=1 # Skip first entry in table
                if self.rename_table.cellWidget(i,1).text() != "": # If rename entry is not empty, save under the entered name
                    plot_name = self.rename_table.cellWidget(i,1).text()
                    canvas.axes.figure.savefig(os.path.join(output_dir,plot_name), bbox_inches = 'tight', dpi=500)
                    pb.setValue(pb.value()+increment) # Update progress bar
                    QtCore.QCoreApplication.processEvents()
                else: # If rename entry is empty, save under original name
                    plot_name = self.rename_table.item(i,0)
                    canvas.axes.figure.savefig(os.path.join(output_dir,plot_name), bbox_inches = 'tight', dpi=500)
                    pb.setValue(pb.value()+increment) # Update progress bar
                    QtCore.QCoreApplication.processEvents()
        self.hide() # Close naming window after done saving figures
        self.label.setVisible(False)
        pb.setValue(0)


    def go_back(self):
        self.hide()
        self.label.setVisible(False)
        self.rename_table.clear()
        


def main():
    if hasattr(sys, 'frozen'):
        global_settings.appdir = sys.executable
        if platform.system() == 'Windows':
            global_settings.appdir = global_settings.appdir[:global_settings.appdir.rfind("\\") + 1]
        else:
            global_settings.appdir = global_settings.appdir[:global_settings.appdir.rfind/("/") + 1]
    else:
        global_settings.appdir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            global_settings.appdir += '\\'
        else:
            global_settings.appdir += '/'

    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = Qt.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("KmCalc")
    global_settings.main_window = MainWindow(os.getcwd())
    global_settings.loading_window = LoadingWindow()
    global_settings.results_window = ResultsWindow()
    print("Current app directory is: " + global_settings.appdir + "\n")
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
