# =========================================================================================
# KmCalc:
# Program for high throughput enzyme kinetic assay data analysis using Michaelis–Menten kinetic model. 
# =========================================================================================




### Imports
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.optimize import curve_fit
from uncertainties import ufloat
import sys
import os, platform
import global_settings as gs
from PyQt5 import QtWidgets, Qt, QtCore, uic
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar

class MainWindow(QtWidgets.QMainWindow):
    """
    Provides a window where the user can set the analysis parameters such as data file, units, etc.

    This class provides a platform where the user can set the analysis parameters. It then performs calculations and analysis
    on the input data to find the reaction rate of each data curve.

    Inputs:
    - Input data file
    - Save output data option (optional)
    - Output file name (optional)
    - Concentration and time units
    - k_cat and K_M initial guesses

    Ouput:
    - Output file (optional)
    - Results graph display
    """
    def __init__(self,appdir):
        """
        Initializes main window and certain parameters.

        Run when:
        - Program is run.
        """
        super(MainWindow, self).__init__()
        uic.loadUi(gs.uidir + "analysis.ui", self)                      # Load main window ui file
        self.mwfg = self.frameGeometry()                                   # Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  # Center window
        self.mwfg.moveCenter(self.cp)                                      # Center window
        self.move(self.mwfg.topLeft())                                     # Center window

        ### Variable Initialization
        self.input_path = ""   # Initialize input file path variable
        self.output_path = ""  # Initialize output file path variable
        self.kcat_guess = 13.7 # Initialize kcat inital guess variable (s^-1)
        self.km_guess = 0.130  # Initialize Km inital guess variable (mM)
        # Initial guesses found here: http://book.bionumbers.org/how-many-reactions-do-enzymes-carry-out-each-second/

        ### Style Modifications
        groupbox_style = """
        QGroupBox:title{subcontrol-origin: margin;
                        left: 10px;
                        padding: 0 5px 0 5px;}
        QGroupBox#Step1{border: 2px solid rgb(68,114,196);
                        border-radius: 5px;
                        font: 14pt "Arial";
                        font: bold;
                        margin-top: 10px;}"""

        self.Step1.setStyleSheet(groupbox_style)
        self.Step2.setStyleSheet(groupbox_style.replace("Step1", "Step2").replace("rgb(68,114,196)","rgb(68,114,196)"))

        ### Button Modifications

        self.progressBar.hide()   # Hide the the progress bar until it is needed
        self.loading_label.hide() # Hide the the loading label until it is needed

        conc_units_list = ['M','mM','uM','nM','pM'] # List of concentration units
        time_units_list = ['s','m','h']             # List of time units
        rate_units_list = []                        # Initialize list for rate units to be made

        for c in conc_units_list:                  # for each concentration unit
            for t in time_units_list:              # for each time unit
                rate_units_list.append(f'{c}/{t}') # Add the new rate unit
        
        self.time_conv_dict = { # Dictionary of convertions from input-time-units to output-time-units (i.e., min * (60 sec / 1 min) = sec)
            's': {'s': 1,'m': 1/60,'h': 1/3600},
            'm': {'s': 60,'m': 1,'h': 1/60},
            'h': {'s': 3600,'m': 60,'h': 1}}
        
        self.conc_conv_dict = {
            'M': {'M': 1,'mM': 10**3,'uM': 10**6,'nM': 10**9,'pM': 10**12},
            'mM': {'M': 10**-3,'mM': 1,'uM': 10**3,'nM': 10**6,'pM': 10**9},
            'uM': {'M': 10**-6,'mM': 10**-3,'uM': 1,'nM': 10**3,'pM': 10**6},
            'nM': {'M': 10**-9,'mM': 10**-6,'uM': 10**-3,'nM': 1,'pM': 10**3},
            'pM': {'M': 10**-12,'mM': 10**-9,'uM': 10**-6,'nM': 10**-3,'pM': 1}}

        self.rate_units_box.addItems(rate_units_list) # Add rate units to rate units selection box
        self.rate_units_box.setCurrentText('mM/s')    # Set defult units to mM/s

        self.kcat_guess_box.setValue(self.kcat_guess) # Set defult kcat guess
        self.km_guess_box.setValue(self.km_guess)     # Set defult Km guess
        self.rate_units_box.currentTextChanged.connect(self.update_guesses) # If the user changes the disired rate units update the displayed initial guesses
        self.kcat_guess_box.valueChanged.connect(lambda: self.kcat_guess == self.kcat_guess_box.value()) # If the kcat guess is changed record the change
        self.km_guess_box.valueChanged.connect(lambda: self.km_guess == self.km_guess_box.value())       # If the Km guess is changed record the change

        self.browse_button.clicked.connect(self.browse_file)         # If browse button is clicked run browse file function
        self.save_data_box.stateChanged.connect(self.disable_output) # If the save data checkbox is clicked check to run appropriate changes
        self.submit_button.clicked.connect(self.submit)              # If the submit button is clicked run the submit fucntion

        # self.show() # Show the window
        self.hide()

    def gather_args(self):
        """
        Collects input parameters.

        Inputs:
        - Input data file path
        - Output data file name (optional)
        - kcat guess
        - Km guess

        Run when:
        - The submit button is pressed.
        """

        if self.output_edit.text() == "": # if no output path is selected...
            self.output_path = os.path.basename(self.input_path).strip('.xlsx') + "-output.xlsx" # Get the input file name, remove the .xlsx, and add "-output" on the end
        elif (self.output_edit.text().endswith(".xlsx") or self.output_edit.text().endswith(".csv")): # if the the user is saving the output as an excel or csv file
            self.output_path = gs.appdir + self.output_edit.text() # Set the output file path as the app directory
        else:
            self.output_path = gs.appdir + self.output_edit.text() + ".xlsx" # Otherise, save the output file in the app directory as an excel file

        self.kcat_guess = self.kcat_guess_box.value() # Get the kcat guess
        self.km_guess = self.km_guess_box.value() # Get the Km guess

    def submit(self):
        """
        Executes analysis, hides main window, and shows results window.

        Run when:
        - Submit button is pressed
        """

        if self.input_path != '': # if the input path is not empty (a file is selected)
            self.gather_args()
            self.find_parameters(self.input_path, self.output_path)
        else:
            self.input_edit.setText('Please select a file to analyze')
    
    def disable_output(self):
        """
        Checks if the save data checkbox is checked and, if so, enables the output file name line edit.

        Run when:
        - The save data checkbox is clicked.
        """
        if self.save_data_box.isChecked():
            self.output_name_label.setDisabled(False)
            self.output_edit.setDisabled(False)
        else:
            self.output_name_label.setDisabled(True)
            self.output_edit.setDisabled(True)

    def update_guesses(self):
        """
        Updates the displayed kcat and Km guesses.

        Run when:
        - Rate units are changed
        """
        self.get_disired_units()
        current_kcat_units = self.kcat_guess_units_label.text().split('/')[-1]
        current_km_units = self.km_guess_units_label.text()
        self.kcat_guess = self.kcat_guess_box.value() / self.time_conv_dict[current_kcat_units][self.disired_time_units]
        self.km_guess = self.km_guess_box.value() * self.conc_conv_dict[current_km_units][self.disired_conc_units]
        self.kcat_guess_box.setValue(self.kcat_guess)
        self.km_guess_box.setValue(self.km_guess)
        self.kcat_guess_units_label.setText(f'1/{self.disired_time_units}')
        self.km_guess_units_label.setText(self.disired_conc_units)

    def browse_file(self):
        """
        Opens file selection dialog.

        Run when:
        - Browse button pressed
        """
        myFile = QtWidgets.QFileDialog.getOpenFileName(caption="Please select an input file...", # Sets the title of the window
                                                        directory=os.getcwd(),                   # Sets the starting directory as the current working directory
                                                        filter='excel (*.xlsx)')                 # Sets the file filter
        
        if (myFile[0] != ""): # if the absolute file path is not empty (a file is selected)
            my_file = str(myFile[0]) # Extracted the file path in a string
            self.input_edit.setText(my_file) # Display the absolute file path
            self.input_path = my_file

            if platform.system() == "Windows": # if the user is using a windows OS
                self.input_path = self.input_path.replace("/","\\") # Use \ instead of /
            else:
                self.input_path = self.input_path.replace("\\","/")

    """ THIS SECTION HOLDS THE FUNCTIONS THAT CALCULATE KINETIC PARAMETERS """

    def find_parameters(self,input_path, output_path):
        """
        Performs analysis on input data.

        Core function of the program.

        Run when:
        - Submit button is pressed
        """

        self.progressBar.show()
        self.loading_label.show()
        self.progressBar.setValue(0) # Set progress bar
        try:
            raw_df, col_map, conversions = self.read_input(input_path) # Get data in dataframe
            self.progressBar.setValue(5) # Set progress bar
            QtCore.QCoreApplication.processEvents() # Update progress Bar
            self.loading_label.setText('Data successfully read')

            ### Find units
            self.get_units(raw_df, col_map, conversions)
            self.progressBar.setValue(20)
            QtCore.QCoreApplication.processEvents()
            self.loading_label.setText('Units identified')

            df = raw_df[[f'Time ({self.time_units})']+list(col_map['Well'])] # Only keep columns noted in column map
            rates_abs_tin, rates_r2, linear_region = self.find_rates(df, col_map) # Get absorbance rates, rate R^2 values, and linear regions
            
            ### Convert to disired units
            rates_cin_tin = rates_abs_tin * conversions.loc['absorbance2concentration', 'value'] # In input-conc-units/input-time-units
            rates_cin = rates_cin_tin / self.time_conv_dict[self.standard_time_units][self.disired_time_units] # In input-conc-units/output-time-units
            rates = rates_cin * self.conc_conv_dict[self.data_conc_units][self.disired_conc_units]

            new_enz_conc_units_values = col_map.loc[:,f'[Enzyme] ({self.enz_conc_units})'] * self.conc_conv_dict[self.enz_conc_units][self.disired_conc_units]
            new_sub_conc_units_values = col_map.loc[:,f'[Substrate] ({self.sub_conc_units})'] * self.conc_conv_dict[self.sub_conc_units][self.disired_conc_units]
            necu_df = pd.DataFrame(new_enz_conc_units_values.values, columns=[f'[Enzyme] ({self.disired_conc_units})'])
            nscu_df = pd.DataFrame(new_sub_conc_units_values.values, columns=[f'[Substrate] ({self.disired_conc_units})'])

            col_map[f'[Enzyme] ({self.enz_conc_units})'] = necu_df[f'[Enzyme] ({self.disired_conc_units})']
            col_map[f'[Substrate] ({self.sub_conc_units})'] = nscu_df[f'[Substrate] ({self.disired_conc_units})']

            col_map_columns_list = list(col_map.columns)
            col_map_columns_list[3:5] = [f'[Substrate] ({self.disired_conc_units})',f'[Enzyme] ({self.disired_conc_units})']
            col_map.columns = col_map_columns_list

            self.progressBar.setValue(60)
            QtCore.QCoreApplication.processEvents()
            self.loading_label.setText('Rates calculated')

            replicates = self.split_replicates_and_subtract_blank(rates, col_map) # Get dictionary of dataframes for replicates for each sample ID ({sample_id: {rep1: conc_rate_df, rep2:...}})
            out_df, self.out_dict = self.get_output_df(replicates) # Get output dataframe and output parameter values (e.g. km, vmax, etc)
            self.rates_out = self.write_report(output_path, out_df, replicates, rates_r2, col_map, linear_region) # Get final dataframe
            self.progressBar.setValue(80)
            QtCore.QCoreApplication.processEvents()

            gs.results_window.generate_plots(self.rates_out, self.out_dict) # Generate plots
            gs.results_window.was_checked()
            self.progressBar.setValue(90)
            QtCore.QCoreApplication.processEvents()
            self.loading_label.setText('Plots generated')
            
            self.progressBar.setValue(100) 
            QtCore.QCoreApplication.processEvents()

            self.progressBar.hide()
            self.loading_label.hide()

            self.hide() # Hide main window
            gs.results_window.show()
        except Exception as e:
            exc_type, error_msg, exc_tb = sys.exc_info() # Extract error information
            print(error_msg) # Print the error message
            self.loading_label.setText(f'ERROR - Line {exc_tb.tb_lineno} - Type {exc_type} - {error_msg}') # Display the error information

    def get_disired_units(self):
        """
        Gets rate, concentration, and time units from rate units selection box.

        Run when:
        - Rate units are changed
        - Submit button pressed
        """
        
        self.rate_units = self.rate_units_box.currentText()
        self.disired_conc_units, self.disired_time_units = self.rate_units.split('/')

    def get_units(self, raw_df, col_map, conversions):
        """
        Extracts units from the data.

        Inputs:
        - Raw data dataframe
        - Data column map dataframe
        - Conversions dataframe

        Outputs:
        - Concentration and times units of input data

        Run when:
        - Submit button pressed
        """
        self.get_disired_units()
        t = raw_df.columns[0] # Gets the column heading for time
        self.time_units = t[t.find('(') + 1: t.find(')')]# Finds the time units

        ### The conc units should be the same but I included them
        # so that potential error checking could be implimented in the future.
        # One variable (e.g. conc_units) could be used instead 
        s_col = col_map.columns[3] # Gets the column heading for [substrate]
        e_col = col_map.columns[4] # Gets the column heading for [enzyme]
        self.sub_conc_units = s_col[s_col.find('(') + 1: s_col.find(')')] # Extracts the substrate concentration units
        self.enz_conc_units = e_col[e_col.find('(') + 1: e_col.find(')')] # Extracts the enzyme concentration units

        ### Establish standard readable time units for later conversions
        if 's' and not('minutes' or 'hours') in self.time_units.lower():
            self.standard_time_units = 's'
        elif 'm' in self.time_units.lower():
            self.standard_time_units = 'm'
        elif 'h' in self.time_units.lower():
            self.standard_time_units = 'h'

        self.data_conc_units, abs1cm_unit = conversions.loc['absorbance2concentration', 'units'].strip().split('/')
    
    def read_input(self,input_path):
        """
        Reads the input file and holds the data in dataframes.

        Input:
        - Input file path

        Output:
        - Raw data dataframe
        - Data column map dataframe
        - Conversions dataframe

        Run when:
        - Submit button is pressed
        """
        raw_df = pd.read_excel(input_path, sheet_name='Raw Data') 
        col_map = pd.read_excel(input_path, sheet_name='Column Map')
        conversions = pd.read_excel(input_path, sheet_name='Conversion Factors', index_col='name')
        return raw_df, col_map, conversions

    def write_report(self, output_path, out_df, replicates, rates_r2, col_map, linear_region):
        """
        Organizes rate values and statistics to be used in results graph and/or to be put in the output file.

        Inputs:
        - Output file path
        - Output dataframe
        - Replicates dictionary
        - Rate R^2 values
        - Column map dataframe
        - Data linear regions

        Output:
        - Final rates dataframe
        - Results output file (optional)

        Run when:
        - Submit button is pressed
        """

        frames = []
        for k,v in replicates.items(): # for each of the replicates dictionaries...
            for k1,v1 in v.items(): # for the reaction rate values in those dictionaries...
                frames.append(v1) # Make a list of the rate values (nested dataframe in a list)
        rates_r = pd.concat(frames)['rate'] # Concatenate the values into a dataframe and extract the rate values
        rates_r2_r = rates_r2
        rates_r.name = f'Rate ({self.rate_units})' 
        rates_out = col_map.set_index('Well').join(rates_r).join(rates_r2_r).join(linear_region) # On the column map dataframe, set the index to the wells (A1, A2, etc), add on the rates, the R^2 values, and the linears regions
        
        if self.save_data_box.isChecked(): # if the user wants to save the output data (i.e. the "Save Data" check box is checked)
            if output_path.endswith('xlsx'): # if the output is an excel file
                with pd.ExcelWriter(output_path, engine='xlsxwriter') as w: # Write to the output excel...
                    out_df.to_excel(w, sheet_name='parameters') # The parameters
                    rates_out.to_excel(w, sheet_name='rates') # The rates and accociated values
                    w.save() # Save the file
                    print(f'Output written to: {os.path.abspath(output_path)}')
            else: # else do the same but for a csv file
                    out_df.to_csv(os.path.join(output_path,'parameters.csv'))
                    rates_out.to_csv(os.path.join(output_path,'rates.csv'))

        return rates_out

    def find_rates(self, df, col_map):
        """
        Calculates the rate from each data curve.

        Uses linear regression of the data in the linear region to get the rates.

        Inputs:
        - Raw data dataframe
        - Column map dataframe

        Outputs:
        - Rates series
        - Rates R^2 value series
        - Linear regions series

        Run when:
        - Submit button is pressed
        """

        time = df[f'Time ({self.time_units})'] # Extract time values
        rates = pd.Series(index=df.columns[1:], name='rate') # Initialize series for rxn rates
        rates_r2 = pd.Series(index=df.columns[1:], name='R2') # Initialize series for rate R^2 values
        linear_region = pd.Series(index=df.columns[1:], name=f'Linear Region ({self.time_units})') # Initialize series for linear regions 

        for col_id in df: # for each column in the dataframe...
            col = df[col_id] # Extract column values
            if col.name == f'Time ({self.time_units})': # if the time column is detected
                continue # skip to the next loop
            
            in_linear_region = self.get_linear_region(time, col, col_map) # get the linear region (boolean array used to slice out the time cooresponding to the linear region in the data)
            result = linregress(x=time[in_linear_region], y=col[in_linear_region]) # Perform a linear regression where the x values are the times in the linear region and the y is the data in the linear region
            rates[col_id] = abs(result[0]) # Extract the rates from the lin. reg. (Note: assumes depletion curves can find product production rate)
            rates_r2[col_id] = result[2] # Extract the R^2 values
            linear_region[col_id] = '{}-{}'.format(time[in_linear_region].iloc[0],time[in_linear_region].iloc[-1]) # Get the start and end time of the linear region and formate it as 'start-end'
        return rates, rates_r2, linear_region

    def get_linear_region(self, time, col, col_map):
        """
        Extracts the linear regions from the column map.

        Inputs:
        - Time values
        - Experimental data values
        - Column map dataframe

        Output:
        - Boolean array with true values in the linear region

        Run when:
        - Submit button is pressed
        """
        start = col_map.loc[ col_map['Well'] == col.name, 'Start Linear Region'].values[0]
        end = col_map.loc[ col_map['Well'] == col.name, 'End Linear Region'].values[0]
        return (time >= start) & (time <= end) # return the boolean array of the union between the time values greater than or equal to the start time and less than or equal to the end time
    
    def split_replicates_and_subtract_blank(self, rates, col_map):
        """
        Organizes the results by sample ID and replicate and subtracts the control rates (samples with no substrate) from the other rate values.

        Inputs:
        - Final rate values
        - Column map dataframe

        Output:
        - Replicates dictionary

        Run when:
        - Submit button is pressed
        """
        # format replicates as: # {sample_id: {rep1: conc_rate_df, rep2:...}}
        replicates = {} # Initialize replicates dictionary
        for sample_id in col_map['Sample ID'].unique(): # for each unique sample ID...
            enz_map = col_map[col_map['Sample ID'] == sample_id] # Get info around the sample
            rep_dict = {} 
            for rep in enz_map['Replicate'].unique(): # for each unique replicate of the sample...
                rep_map = enz_map[enz_map['Replicate'] == rep] # Get the info around the replicate
                rep_df = rep_map.set_index('Well').join(rates) # Add rates
                # subtract blank
                blank_value = rep_df[rep_df[f'[Substrate] ({self.disired_conc_units})'] == 0]['rate'].values # Find the rate for the well with no substrate
                rep_df['rate'] = rep_df['rate'] - blank_value # Substract the blank well rate to get the pure rates
                rep_dict[rep] = rep_df # Set the info for each replicate in the replicate dictionary
            replicates[sample_id] = rep_dict # Set the replicate dictionary for each sample ID
        return replicates 

    def get_output_df(self,replicates):
        """
        Creates a dictionary and dataframe with the enzyme kinetics parameters for each sample ID.

        Inputs:
        - Replicates dictionary

        Outputs:
        - Output parameters dictionary
        - Output parameters dataframe

        Run when:
        - Submit button is pressed
        """
        out_dict = {}
        for sample_id, reps in replicates.items():
            # From a statistical point of view it makes more sense to fit the model to all available data instead of
            # treating replicates individually, thus:
            all_rep_df = pd.concat([rep_df for id, rep_df in reps.items()]) # Concatenate all of the data for the replicates of a sample together
            km, vmax = self.find_km_vmax(all_rep_df) # Calculate km and vmax
            # Enzyme conc (ug/mL)
            kcat = vmax/all_rep_df.loc[ # Calculate kcat (vmax/[total enzyme])
                all_rep_df.index[0], # For the first well in the sample id
                f'[Enzyme] ({self.disired_conc_units})'] # What is the [enzyme] Note: assumes [enzyme] is constant across the wells for a sample 
            kcatkm = kcat/km # Calculate enzymatic efficiency
            out_dict[sample_id] = {f'km ({self.disired_conc_units})':km,
                                    f'vmax ({self.rate_units})':vmax,
                                    f'kcat (1/{self.disired_time_units})':kcat,
                                    f'kcat/km (1/{self.disired_time_units}/{self.disired_conc_units})':kcatkm} # Set dictionary with values
        out_df = pd.DataFrame(out_dict).transpose()

        return out_df, out_dict

    def mm_model(self,s, vmax, km):
        """
        Michaelis–Menten equation.

        Inputs:
        - Substrate concentration
        - Vmax
        - Km

        Output:
        - Rate of reaction

        Run when:
        - Submit button pressed (for curve fitting)
        - Plotting Fitted curve
        """
        return (vmax*s)/(km + s)

    def find_km_vmax(self,rep_df):
        """
        Uses robust non-linear regression to get kinetic parameters.

        The paremeter covariance estimate procedure is explained here: https://stackoverflow.com/questions/14854339/in-scipy-how-and-why-does-curve-fit-calculate-the-covariance-of-the-parameter-es

        Bounds and inital guesses on km and vmax parameters (concentration units), from Cell Biology by the Numbers: http://book.bionumbers.org/how-many-reactions-do-enzymes-carry-out-each-second/
        
        Input:
        - Replicate dataframe

        Output:
        - Km value
        - Vmax value

        Run when:
        - Submit button is pressed
        """

        mean_Etot = rep_df[f'[Enzyme] ({self.disired_conc_units})'].mean() # Get total enzyme concentration (should all be the same)
        bounds = (np.full([2,], 0), np.array([mean_Etot * (10e7) / self.time_conv_dict['s'][self.disired_time_units] , 10e7 * self.conc_conv_dict['uM'][self.disired_conc_units]])) # Makes a tuple with 2 arrays of size 2 filled with 0 and 10e4 respectively # WARNING: THIS BOUNDS ASSUME CERTAIN INPUT UNITS! # @ change so assumption doesn't need to be made
        x0 = np.array([mean_Etot * self.kcat_guess / self.time_conv_dict['s'][self.disired_time_units], self.km_guess * self.conc_conv_dict['uM'][self.disired_conc_units]]) # Initial guesses for Vmax and Km respectively
        
        # @ Changing units changes km value for some reason
        popt, pcov = curve_fit(self.mm_model, # Equation to be fit
                            xdata=rep_df[f'[Substrate] ({self.disired_conc_units})'], # x-values ([substrate])
                            ydata=rep_df['rate'], # y-values (rate of reaction)
                            p0=x0, # Set initial guesses for parameters
                            bounds=bounds, # tuple of arrays defining the bounds of the values as being between 0 and 10e4
                            loss='soft_l1') # From documentation, "‘soft_l1’ : rho(z) = 2 * ((1 + z)**0.5 - 1). The smooth approximation of l1 (absolute value) loss. Usually a good choice for robust least squares."
        perr = np.sqrt(np.diag(pcov)) # Standard deviation of the parameters. See curve_fit documentation
        vmax = ufloat(popt[0], perr[0]) # Saves the value with +/- uncertainty i.e popt[0] +/- perr[0]
        km = ufloat(popt[1], perr[1])

        return km, vmax

class ResultsWindow(QtWidgets.QMainWindow):
    """
    Generates results plot.

    Creates a graph that displays data points with error bars, fit Michaelis–Menten equation curve, and legend with kinetic parameter values (kcat, Km, Vmax) for each sample ID.

    Input:
    - Analysis results from main window

    Output:
    - Results display
    """
    def __init__(self):
        """
        Initializes results window and certain parameters.

        Run when:
        - Program is run
        """
        super(ResultsWindow, self).__init__()
        uic.loadUi(gs.uidir + "results_window.ui", self)
        self.plots_generated = False
        self.mwfg = self.frameGeometry()  # Center window
        self.cp = QtWidgets.QDesktopWidget().availableGeometry().center()  # Center window
        self.mwfg.moveCenter(self.cp)  # Center window
        self.move(self.mwfg.topLeft())  # Center window
        self.window_width = self.frameGeometry().width()

        ### Initialize layouts for figures
        self.widget_layout = QtWidgets.QVBoxLayout()
        self.widget_layout.setContentsMargins(0,0,0,0)
        self.widget_2_layout = QtWidgets.QVBoxLayout()
        self.widget_2_layout.setContentsMargins(0,0,0,0)

        ### Connect functions to buttons
        self.back_button.clicked.connect(self.go_back)

        self.checkboxes_made = False

        groupbox_style = """
                QGroupBox:title{subcontrol-origin: margin;
                                left: 10px;
                                padding: 0 5px 0 5px;}
                QGroupBox#all_box{border: 2px solid rgb(68,114,196);
                                border-radius: 5px;
                                font: 14pt "Arial";
                                font: bold;
                                margin-top: 10px;}"""

        self.all_box.setStyleSheet(groupbox_style)

    def was_checked(self):
        """
        Establishes a connection between the replot function and the state of the sample ID checkboxes.

        Run when:
        - Submit button is pressed
        """
        for sample, box in self.checkbox_dict.items(): # for each sample ID checkbox
            box.stateChanged.connect(self.replot) # Establish a connection to the replot function

    def replot(self):
        """
        Clears the plot of the previous graph and calls the generate plots function to replot with the sample IDs that are checked.

        Run when:
        - A sample ID checkbox is clicked
        """
        self.plots_generated = False
        self.widget_canvas.deleteLater() # Delete the previous canvas widget
        self.toolbar.deleteLater() # Delete the previous toolbar widget
        self.generate_plots(gs.main_window.rates_out,gs.main_window.out_dict)

    def go_back(self, return_home = False):
        """
        Clears the results window of previous widgets and layouts along with reverting certain indicator variables.

        Prepares the results window to display a different set of results.

        Run when:
        - Back button is pressed
        - Return home button is pressed
        """
        self.hide() # Hide the results window
        self.plots_generated = False # Reset plots generated indicator variable
        self.checkboxes_made = False # Reset checkboxes made indicator variable
        self.widget_canvas.deleteLater() # Delete the previous canvas widget
        self.toolbar.deleteLater() # Delete the previous toolbar widget
        self.scroll_layout.deleteLater() # Delete the scroll area layout
        for id, checkbox in self.checkbox_dict.items():
            checkbox.deleteLater() # Delete each checkbox
        self.checkbox_dict = {} # Reset checkbox dictionary
        
        if not return_home: # if the user is not returning to the home window
            gs.main_window.show() # Show the main window again
    
    def generate_plots(self, rates_out, out_dict):
        """
        Generates the plot and checkboxes, if not already made, on the results window.

        This is the main function of the results window.

        Run when:
        - Submit button is clicked
        - A sample ID checkbox is clicked
        """
        # Note: The lists below may be useful later but are currently unused
        self.sample_ids = []
        self.enz_ids = []
        self.sub_ids = []

        ### Clear out old widgets in layout
        self.widget_canvas = MplCanvas(self, width=5, height=4, dpi=100) # Initialize new Canvas
        self.toolbar = NavigationToolbar(self.widget_canvas, self) # Initialize toolbar for multisample plot
        self.widget_layout.addWidget(self.toolbar) # Add canvas to widget 1 layout
        self.widget_layout.addWidget(self.widget_canvas) # Add canvas to widget 1 layout
        self.multiwidget.setLayout(self.widget_layout) # Add canvas (with layout) to widget 1 object

        if not self.checkboxes_made: # if the checkboxes haven't already been made
            self.scroll_layout = QtWidgets.QGridLayout() # Establish scroll area layout
            self.scroll_area.setLayout(self.scroll_layout) # Add the layout to the scroll area
            self.scroll_area.adjustSize() # Adjust the size of the checkbox to fit the text
            self.checkbox_dict = {} # Initialize dictionary to keep track of checkboxes
            for i, contents in enumerate(out_dict.items()): # for each sample ID...
                sample_checkbox = QtWidgets.QCheckBox(contents[0].replace('~',' ')) # Create a checkbox for a sample ID 
                sample_checkbox.adjustSize() # Adjust the size of the checkbox to fit the text
                sample_checkbox.setChecked(True) # Set the checkbox to checked
                self.checkbox_dict[contents[0].replace('~',' ')] = sample_checkbox # Save the checkbox in the dictionary
                self.scroll_layout.addWidget(sample_checkbox,i,0) # Add the checkbox to the scroll area
                
                checkbox_width = sample_checkbox.geometry().width() # Get the width of the checkbox
                if checkbox_width >= self.scroll_area.geometry().width(): # if the width of the checkbox is greater than the width of the scroll area
                    self.scroll_area.setMinimumWidth(int(checkbox_width+30)) # Set the minimum width of the scroll area to the width of the checkbox with some buffer                    

            self.checkboxes_made = True # Checkboxes have been made


        counter = -1 # Counter for selecting plot parameters (color, point shapes, etc.)
        plot_styles = ["o","v","^","s","D","P","8","H","X","*","p"]
        colors = ["cornflowerblue","indianred","gold","seagreen","rebeccapurple","orangered","lightslategrey","cadetblue","steelblue","navy","dimgray"]
        for sample_id, params in out_dict.items(): # for each sample ID and its parameters
            counter +=1 # Iterate counter. Put above 'continue' command to retain data colors on graph
            self.sample_ids.append(sample_id) # Add sample ID to list

            if not self.checkbox_dict[sample_id.replace('~',' ')].isChecked(): # if the sample ID checkbox is not checked
                continue # Skip to the next sample ID

            sepi = sample_id.find('~') # Separator index: finds the location of the ~ separator
            self.enz_ids.append(sample_id[0:sepi]) # Get the enzyme ID
            self.sub_ids.append(sample_id[sepi+1:len(sample_id)]) # Get the substrate ID

            erates = rates_out[rates_out['Sample ID'] == sample_id] # Extracts reaction rates for the sample ID
            
            ### Generate figures
            subs = erates[f'[Substrate] ({gs.main_window.disired_conc_units})'].unique() # Get each unique substrate concentration

            for i, sub in enumerate(subs): # for each unique substrate...
                rerates = erates[erates[f'[Substrate] ({gs.main_window.disired_conc_units})']==sub] # Get the information for the substrate concentration
                # Get plot variables
                rates = rerates[f'Rate ({gs.main_window.rate_units})'] # Get rates from the substrate concentration
                e_rates = rates.mean() # Get rates average
                s_conc = rerates[f'[Substrate] ({gs.main_window.disired_conc_units})'].mean() # Get concentrations (taking mean just conveniently puts same concentration in 1 value)
                rerror = rates.std() # Find the standard deviation of the rates
                self.widget_canvas.axes.errorbar(s_conc, e_rates,yerr = rerror, ecolor = colors[counter], capsize = 3, marker=plot_styles[counter],c=colors[counter]) # Add averaged data points to total plot with std error bars

            self.ps_conc = np.linspace(min(subs), max(subs)*1.2, 300) # Create x-values for fitted curve plot
            km = params[f'km ({gs.main_window.disired_conc_units})'] # Get the Km value
            vmax = params[f'vmax ({gs.main_window.rate_units})'] # Get the Vmax value
            kcat = params[f'kcat (1/{gs.main_window.disired_time_units})'] # Get the kcat value
            self.pe_rates = [gs.main_window.mm_model(x, km=km.nominal_value, vmax=vmax.nominal_value) for x in self.ps_conc] # Create the y-values for the fitted curve plot

            ### Plot enzyme's data on total figure
            self.widget_canvas.axes.plot(self.ps_conc, self.pe_rates, c=colors[counter], ls='--',
            label=(r'$\bf{}:$'.format("{" + sample_id.replace("_","-") + "}")+'\nK$_m$ ({}) = {:.3fP}\nV$_{{max}}$ ({}) = {:.3fP}\nk$_{{cat}}$ (1/{}) = {:.3fP}'.format(gs.main_window.disired_conc_units,km,gs.main_window.rate_units,vmax,gs.main_window.disired_time_units,kcat))) # Create plot and add kinetic parameters to the legend


            ### Total plot details
            self.widget_canvas.axes.legend(fontsize=8).set_draggable(state=True)
            self.widget_canvas.axes.set_xlabel(f'Initial substrate ({gs.main_window.disired_conc_units})')
            self.widget_canvas.axes.set_ylabel(f'Rate ({gs.main_window.rate_units})')

            self.widget_canvas.axes.set_ylim([0, max(rates_out[f'Rate ({gs.main_window.rate_units})'])*1.1]) # Sets the y-axes limit to slightly above the largest rate value
            self.widget_canvas.axes.set_xlim([0, max(subs*1.1)]) # Sets the x-axes limit to slightly above the largest substrate value



        self.plots_generated = True # Set bool for whether plots have been generated to true.

class MplCanvas(FigureCanvasQTAgg):
    """
    Used to initialize canvas figure for plotting.
    """
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi,tight_layout=True) # Initialize instance of a figure
        self.axes = fig.add_subplot(111) # Add a subplot to the figure with 1 row, 1 column, and indexed at the first position (i.e., just add a plot)
        self.axes.clear() # Clears the figure
        super(MplCanvas, self).__init__(fig)



def main():
    """
    Executes program.
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
        gs.appdir = os.path.dirname(os.path.abspath(__file__)) # appdir used to keep track of the app directory
        gs.uidir = os.path.dirname(os.path.abspath(__file__)) # uidir used to keep track of the ui folder directory
        if platform.system() == 'Windows':
            gs.appdir += '\\'
            gs.uidir += "\\ui_files\\"
        else:
            gs.appdir += '/'
            gs.uidir += "/ui_files/"


    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = Qt.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("KmCalc-analysis")

    if platform.system() == 'Windows':
        gs.slashreplace = '\\'
    else:
        gs.slashreplace = '/'

    gs.main_window = MainWindow(os.getcwd())
    gs.results_window = ResultsWindow()
    gs.main_window.show()
    print("Current app directory is: " + gs.appdir + "\n")
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()