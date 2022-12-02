"Imports"
### Computer systems
import os, platform
import sys
import global_settings as gs
### Data and number handling
import numpy as np
import pandas as pd
import datetime as dt
### Graphical user interface (PyQt5)
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets, Qt, QtCore, uic

"GUI Classes"
'Invalid row selections pop-up dialog'
class rows_errDialog(QtWidgets.QDialog):
### This dialog is presented when more or less than exactly 2 rows are slected
    def __init__(self):                                         # Initialize class
        super(rows_errDialog, self).__init__()                  # Call the inherited classes __init__ method
        uic.loadUi(gs.uidir + '2rows-dialog.ui', self)                     # Load the .ui file
        self.setStyleSheet('QDialog{background-color: white;}') # Set background color white
        self.ok_btn.clicked.connect(self.close)                 # When 'ok' button is clicked, close the pop-up

'Correct data confirmation pop-up dialog'
class DataConf(QtWidgets.QMainWindow):
### This dialog is presented after the column removal step
# and allows the user to confirm that the data has been
# processed correctly before continuing
    def __init__(self):
        super(DataConf, self).__init__()
        uic.loadUi(gs.uidir + 'data-conf.ui', self)
        self.setStyleSheet('QMainWindow{background-color: white;}')
        self.no_btn.clicked.connect(self.close)

'Raw data cleaning system'
class DataSelection(QtWidgets.QMainWindow):
### This class provides a system for the user to (in order):
# 1. upload an xlsx file of raw data,
# 2. pick out the rows containing the data of interest,
# 3. remove any unnecessary columns.
# Finally, the data is automattically cleaned
# and the DataConf window is prompted for the user
# to confirm that the data has been cleaned correctly
    def __init__(self,appdir):
        super(DataSelection, self).__init__()
        uic.loadUi(gs.uidir + 'data-selection.ui', self)
        self.setStyleSheet('QMainWindow{background-color: white;}')

        self.Bslctc_label.hide() # Initially hides label used in the 2nd step (column removal)
        self.Aprgbar.hide()      # Initially hides progress bar

        self.Aselect_file_btn.clicked.connect(self.selectFileDialog) # Prompts file selection if the 'Select File' button is clicked
        self.submit_btn.clicked.connect(self.submit)                 # Moves to next step if 'Submit' button is clicked
        self.reset_btn.clicked.connect(self.reset)                   # If the user makes a mistake, they can go back to
                                                                        # the original row selection by clicking the 'Reset' button
        gs.colMap.format_btn.clicked.connect(self.format)          # Prompts the user to save a file of the formatted data
        gs.colMap.back_btn.clicked.connect(self.back)              # Goes back to the data selection window and hides the column map window

        # self.show() # Show the window
        self.hide()

    ### Prompts the user  to select a file
    def selectFileDialog(self):
        file_filter = 'excel (*.xlsx)' # File filter has the label 'excel' and only calls .xlsx files
        ### Prompts file selection and saves the file path and file type in a tuple
        response = QFileDialog.getOpenFileNames(caption='Select a data file', # Sets the title of the window
                                                directory=gs.appdir,        # Sets the starting directory as the working directory
                                                filter='excel (*.xlsx)')           # Sets the file filter
        try:
            self.file_path = response[0][0]            # response returns a list inside a tuple => [0][0] indexes the file path out
            self.AlineEdit.setText(self.file_path)     # Sets the file path edit line text to the file path
            self.getRawData()                          # Transfers the contents of the file into a pandas dataframe and upload that data to the table
        except IndexError:
            self.AlineEdit.setText('Please select a file') # If no file is selected the file path edit line displays 'Please select a file'

    ### Brings the user back to the beginning of data selection
    def reset(self):
        self.table.selectionModel().clearSelection() # Set all cells to being not selected

        if not (self.Aslctf_label.isVisible()): # if a label from the row selection step is not visible...
            self.Bslctc_label.hide() # Hide the label in the column selection step...
            ### and show all of the labels in the initial row selection step
            self.Aslctf_label.show()
            self.AlineEdit.show()
            self.Aselect_file_btn.show()
            self.Aslctr_label.show()
            self.table.setSelectionBehavior(QAbstractItemView.SelectRows) # Change the selection mode back to row selection

            self.getRawData() # Transfer the original data onto the table again

    ### Transfers the contents of the file into a pandas dataframe and uploads that data to the table
    def getRawData(self):
        self.df = pd.read_excel(self.file_path, header=None) # Extracts the contents of the excel file and puts them in a dataframe
        self.nrows = self.df.shape[0]                        # Get the number of rows (used in selectRows)
        self.ncolumns = self.df.shape[1]                     # Gets the number of columns (used in selectRows)
        self.dataToTable(self.df, self.table)                # Transfers the contents of the dataframe to the table

    ### Transfers the contents of a given dataframe to a table
    def dataToTable(self, dataframe, table_name):
        num_rows = dataframe.shape[0]                 # Get number of rows
        num_cols = dataframe.shape[1]                 # Get number of columns
        table_name.setRowCount(num_rows)              # Sets the number of rows in the table to the rows measured from the dataframe
        table_name.setColumnCount(num_cols)           # Sets the number of columns in the table to the columns measured from the dataframe
        self.Aprgbar.setRange(0, num_rows * num_cols) # Sets the range of the progress bar from 0 to the number of cells in the table
        self.Aprgbar.show()                           # Show the progress bar

        prgv = 0                            # Initialize progress bar value
        for row in range(num_rows):         # for each row...
            for col in range(num_cols):     # for each column in the row...
                table_name.setItem( row, col, QTableWidgetItem(str(dataframe.iloc[row,col]))) # Set the contents of the cell to
                                                                                                 # the corresponding data in the dataframe
                prgv += 1                   # Increment the p.b. value
                self.Aprgbar.setValue(prgv) # Show the p.b. value on the progress bar
        self.Aprgbar.hide()                 # Once all the data is transfered, hide the progress bar

    ### Gets the 2 user selected rows and removes the rows above and below the data in between (keeping the selected rows)
    def selectRows(self):
        rows = {}
        row_list = self.table.selectedItems() # Creates a list of table widget objects with the information of each selected cell (cells = total columns * number of rows selected)
        try:
            for r in range(int(len(row_list) / self.ncolumns)):  # for the number of rows (number of rows selected = cells / total columns)...
                rows[str(r)] = row_list[self.ncolumns * r].row() # Get the number of the rows selected and hold them in the dictionary as '0' and '1'
            if len(rows) > 2:     # If more that 2 rows are selected
                raise Exception() # Raise an exception (error)
            self.frow = rows['0'] # First row
            self.lrow = rows['1'] # Last row

            ### When a row is removed from a table widget, all of the rows below it move up
            # and become row number r - 1 ( e.g (delete row 1) => row 3 -> row 2).
            # So by repeatedly removing the the same row, one can delete a block of rows.
            for afrow in range(self.frow):                      # for all of the rows above the first selected...
                self.table.removeRow(0)                         # Remove the first row
            for blrow in range(self.nrows-(self.lrow+1)):       # for all of the rows below the last selected...
                self.table.removeRow((self.lrow - self.frow)+1) # Remove the row directly below the last

            self.new_nrows = (self.lrow - self.frow)+1          # Set the new number of rows (rows including and between those selected)
        except:                                          # If more or less than exactly 2 rows selected...
            self.table.selectionModel().clearSelection() # Clear the rows selected (if any)
            gs.rows_errDialog.show()                        # Prompt error pop-up

    ### Gets the unwanted user selected columns and removes them
    def selectColumns(self):
        self.columns = []
        column_list = self.table.selectedItems() # Creates a list of table widget objects with the information of each selected cell

        for c in range(int(len(column_list) / self.new_nrows)):           # for the number of columns selected...
            self.columns.append(column_list[self.new_nrows * c].column()) # Add the number of each selected column to the list

    ### Does final data processing and cleaning
    def finDataProc(self):
        row_range = range(self.table.rowCount())                                 # Creates a range of the number of rows in the table
        col_range = range(self.table.columnCount())                              # Creates a range of the number of columns in the table
        predf = pd.DataFrame({}, index=list(row_range), columns=list(col_range)) # Initialize empty dataframe with the same size as the table

        for col in col_range:                                         # for each column...
            for row in row_range:                                     # for each row...
                predf.iloc[row,col] = self.table.item(row,col).text() # Transfer the corresposnding data from the table to the dataframe

        ### Empty cells transfer to the dataframe as the string object 'nan'.
        # So for them to easily be removed with pandas they need to be
        # converted into numpy not a number (NaN) objects
        nantoNaN = predf.replace('nan',np.nan)
        datadrop = nantoNaN.drop(columns=self.columns)            # Remove the previously selected columns
        self.datafin = datadrop.dropna(axis='columns', how='any') # Remove the columns with non-numbers (NaN)

        ###Convert time to decimal hours format (only works for <24h run) if in h:m:s format.
        # Assumes that the first row contains the data headers.
        # Does not assume that the time data is in the first column and can find the time column
        # as long as it contains the word "time" (e.g time, Time, Time (min), time (s), etc. are acceptable).
        # Simply using "t" as the key word is not used to avoid confusion with potential "TX" (e.g. T2) wells.

        self.headers = self.datafin.iloc[0,:]     # Represent the headers as the contents of the first row
        for loc, item in enumerate(self.headers): # for the location and content of each item in the list...
            self.location = loc                   # Set the location
            if 'time' in item.lower():            # if the item contains "[time]"...
                break                             # Break the loop, setting the location of the time column

        tuv = 0  # Time units variable: indicates if the time units have been identified

        ### Identifies the units of the 'time' column
        if ('hour' or 'hr') in self.headers.iloc[self.location].lower():
            self.timeunits = 'Hours'
            tuv = -1
        if 'min' in self.headers.iloc[self.location].lower():
            self.timeunits = 'Minutes'
            tuv = -1
        if 'sec' in self.headers.iloc[self.location].lower():
            self.timeunits = 'Seconds'
            tuv = -1

        for i, item in enumerate(self.datafin.iloc[:,self.location]): # for the index and content of each item in the time column...
            try:
                if ':' in str(item):
                    temp = dt.time.fromisoformat(str(self.datafin.iloc[i,self.location]))                      # Extract the time information from the item
                    self.datafin.iloc[i,self.location] = float(temp.hour*60 + temp.minute +(temp.second / 60)) # Convert the format into a decimal
                    tuv = 1
            except ValueError:   # If the time data is not in isoformat
                x = 'Do nothing' # Do nothing and skip the above code

        if tuv == 1:
            self.timeunits = 'Minutes'
        elif tuv == 0:
            self.timeunits = 'unknown'

        self.datafin.iloc[0,self.location] = f'Time ({self.timeunits})' # Set a Standardized, potentially updated, time column header
        gs.colMap.initial_timeunits_label.setText(self.timeunits)     # Shows the user the identified time units in the column map selection window
        self.headers = self.datafin.iloc[0, :]
        self.dataToTable(self.datafin, gs.DataConf.conf_table) # Transfer the final data to the confirmation table
        self.datafin.columns = self.headers                      # Set the headers to the contents of the first row
        self.datafin.index = self.datafin.iloc[:,self.location]       # Set the index to the time column
    
    ### Progresses the process to the next step
    def submit(self):
        ### After the 'Sumbit' button is pressed...
        
        ### Row selection step
        if self.Aslctf_label.isVisible(): # if the label associated with the row selection step is visible...
            self.selectRows()             # Run selectRows process

        ### Column selection step
        if self.Bslctc_label.isVisible():  # if the label associated with the column selection step is visible...
            self.selectColumns()           # Run selectColumns process
            self.finDataProc()             # Do final data processing
            gs.DataConf.show()           # Present data confirmation dialog
            gs.DataConf.yes_btn.clicked.connect(
                self.hidetransfer)  # If the 'Yes' button on the data confirmation dialog is clicked,
                                       # transfer to the column map selection process

        if not(gs.rows_errDialog.isVisible()) and not(self.Bslctc_label.isVisible()): # if the row error dialog isn't present
                                                                                      # and the label associated with the
            ### Hide the widgets associated with the row selection                    # column map selection process isn't visible
            self.Aslctf_label.hide()
            self.AlineEdit.hide()
            self.Aselect_file_btn.hide()
            self.Aslctr_label.hide()
            
            self.Bslctc_label.show()                                         # Show the label associated with the row selection process
            self.table.selectionModel().clearSelection()                     # Clear the selected cells
            self.table.setSelectionBehavior(QAbstractItemView.SelectColumns) # Switch the selection mode to select columns

    ### Hides the data selection and confirmation windows and transfers to the column map selection process
    def hidetransfer(self):
        for h, head in enumerate(self.headers):                                  # for the headers of the data...
            gs.colMap.headers_table.setItem(0, h, QTableWidgetItem(str(head))) # Transfer them to the headers table
        gs.DataConf.hide() # Hide the data confirmation dialog
        self.hide()          # Hide the data selection window
        gs.colMap.showMaximized()   # Show the column map selection window

    ### Allows the user to go back to the data selection from the column map selection
    def back(self):
        self.showMaximized()
        gs.colMap.hide()

    ### Formats the raw and column map data into an excel file
    def format(self):
        ### NOTE: This process occurs after the column map selection process. (See colMap class)
        # This process occurs in the DataSelection class for ease of data transfer between windows
        col_range = range(5)        # Sets column range (1-Well 2-Replicate 3-Sample-ID 4-[Sub] 5-[Enz])
        colmapdf = pd.DataFrame({}) # Initailie empty dataframe to store column map data

        for col in col_range: # for each of the columns...
            row = 0           # Starting at the first row...
            while str(gs.colMap.maintable.item(row, col)) != 'None':               # while the cell is not empty...
                colmapdf.loc[row, col] = gs.colMap.maintable.item(row, col).text() # Transfer the contents of the main column map table to a dataframe
                row += 1                                                             # Move to the next row number

        concunits = gs.colMap.conc_combobox.currentText() # Gets concentration units
        totime = gs.colMap.to_combobox.currentText()      # Gets the time units to convert to

        convf = {'Seconds to Minutes': 60**-1,
                 'Seconds to Hours': 3600**-1,
                 'Minutes to Seconds': 60,
                 'Minutes to Hours': 60**-1,
                 'Hours to Seconds': 3600,
                 'Hours to Minutes': 60
        }

        if (self.timeunits != totime) and (self.timeunits != 'unknown'):
            for i, item in enumerate(self.datafin.iloc[:,self.location]):
                try:
                    self.datafin.iloc[i, self.location] = float(self.datafin.iloc[i, self.location]) * convf[f'{self.timeunits} to {totime}']
                except:
                    x = 'Do Nothing'

        self.timeunits = totime
        self.datafin.iloc[0, self.location] = f'Time ({self.timeunits})'

        try:
            ### Prompts file saving
            name = QFileDialog.getSaveFileName(self, caption='Save As...', # Set the window caption
                                               directory=gs.appdir,      # Set the working directory
                                               filter='EXCEL (*.xlsx)'     # Sets the file filter to excel files
                                               )[0]                        # Indexes out the file name
            
            writer = pd.ExcelWriter(name,engine='openpyxl',options={'strings_to_formulas': False}) # Initialize ExcelWriter class for excel file creation and editing

            self.datafin.to_excel(writer,                # Transfer cleaned data to an excel file...
                                  sheet_name='Raw Data', # in a sheet titled 'Raw Data'
                                  index=None,            # Don't include index
                                  header=None)           # Don't include headers

            colmapdf.to_excel(writer,                    # Transfer column map data to the excel file...
                              sheet_name='Column Map',   # in a sheet titled 'Column Map'
                              index=None,                # Don't include index
                              header=['Well', 'Replicate', 'Sample ID', f'[Substrate] ({concunits})', f'[Enzyme] ({concunits})']) # Set the headers

            abtoconcf = gs.colMap.abs_to_conc_spinbox.value() # Get absorbance to concentration conversion factor defined by user
            ### Set conversion factors dataframe (used to convert units in analysis)

            convfactorsdf = pd.DataFrame([['name', 'value', 'units'],
                                          ['absorbance2concentration', f'{abtoconcf}', f'{concunits}/abs*1cm']])

            convfactorsdf.to_excel(writer,                          # Transfer conversion factors data to the excel file...
                                   sheet_name='Conversion Factors', # in a sheet titled 'Conversion Factors'
                                   index=None,
                                   header=None)

            writer.save() # Save the file and edits
        except:
            'Go Back' # Go back and do nothing

'Column Map selection'
class colMap(QtWidgets.QMainWindow):
### The colMap class provide the user with a streamlined way to quickly format the column map.
# This column map acts as a key or legend to the experimental conditions of each well.
# The information includes:
# Well: The assay well (e.g. A1, B5, C8, etc)
# Replicate: The experimental replicate
# Sample ID: Shows the enzyme and substrate used
# [Substrate]: Substrate concentration
# [Enzyme]: Enzyme concentration
    def __init__(self):
        super(colMap, self).__init__()
        uic.loadUi(gs.uidir + 'column-map-selection.ui', self)
        self.setStyleSheet('QMainWindow{background-color: white;}')
        self.savedrownums = [] # Initializes list to track additions or removal of formatted blocks from the set-table
        self.netaddedrows = 0  # Initializes int variable to keep track of net number of rows added or removed
        self.settable.itemChanged.connect(self.tableupdate) # Automatically puts the replicate and sample ID next to any wells added to the set table
        self.repl_spinbox.valueChanged.connect(lambda *args: self.update_repl())      # Updates replicate number if changed
        self.enzid_lineedit.textChanged.connect(lambda *args: self.sampleid_update()) # Updates sample ID if enzyme ID is changed
        self.subid_lineedit.textChanged.connect(lambda *args: self.sampleid_update()) # Updates sample ID if substrate ID is changed
        self.addslct_btn.clicked.connect(self.addwells) # Adds selected wells to the set-table
        self.clear_btn.clicked.connect(self.headers_table.selectionModel().clearSelection) # Clears headers-table selections
        self.cleartab_btn.clicked.connect(self.settable.clearContents) # Clears the contents of the set-table
        self.clearmaintab_btn.clicked.connect(self.clearmaintab)       # Clears the contents of the main table
        ### Clears the selected cells in either set table or main table if the other is clicked
        # so no cell is unintentionally effected by keyboard shortcuts
        self.settable.cellClicked.connect(self.maintable.selectionModel().clearSelection)
        self.maintable.cellClicked.connect(self.settable.selectionModel().clearSelection)
        self.remove_btn.clicked.connect(self.remove)  # Removes most recently added block of data from the set table
        self.add_btn.clicked.connect(self.addtotable) # Appends the contents of in the set table to the end of the contents in the main table
        self.insert_btn.clicked.connect(self.insert)  # Inserts an empty row above the selected row
        self.delete_btn.clicked.connect(self.delete)  # Deletes the selected row

        self.clipboard = QGuiApplication.clipboard() # Initializes class that allows the user to copy and paste
        ### Initaialize keyboard shortcuts
        self.shortcut1 = QShortcut(QKeySequence('delete'), self)
        self.shortcut2 = QShortcut(QKeySequence('backspace'), self)
        self.shortcut3 = QShortcut(QKeySequence('Ctrl+c'), self)
        self.shortcut4 = QShortcut(QKeySequence('Ctrl+v'), self)
        ### Activates keyboard shortcuts
        self.shortcut1.activated.connect(self.delshortcut)
        self.shortcut2.activated.connect(self.delshortcut)
        self.shortcut3.activated.connect(self.copy)
        self.shortcut4.activated.connect(self.paste)
        self.hide()

    ### Pastes text from the clipboard to the selected cells
    def paste(self):
        ### If set table cell(s) are selected
        items_selected = self.settable.selectedItems() # Get list of selected cells
        for item in items_selected: # for each item in the list...
            self.settable.setItem(item.row(), item.column(), QTableWidgetItem(self.clipboard.text())) # Set the contents of the cell to the text from the clipboard
        ### If main table cell(s) are selected...
        items_selected = self.maintable.selectedItems()
        for item in items_selected:
            self.maintable.setItem(item.row(), item.column(), QTableWidgetItem(self.clipboard.text()))

    ### Copies text from a selected cell onto the clipboard
    def copy(self):
        ### If set table cell(s) are selected...
        items_selected = self.settable.selectedItems()                             # Get list of selected cells
        for i in items_selected:                                                   # for each item in the list...
            self.clipboard.setText(self.settable.item(i.row(), i.column()).text()) # Save the cell text to the clipboard

        ### If main table cell(s) are selected...
        items_selected = self.maintable.selectedItems()
        for i in items_selected:
            self.clipboard.setText(self.maintable.item(i.row(), i.column()).text())

    ### Clears the contents of selected cell(s)
    def delshortcut(self):
        items_selected = self.settable.selectedItems() # Get list of selected cell(s)
        ### If a user is deleting a well, replicate, and/or sample ID from the set table
        # they are likely wanting to just delete the whole row.
        # The next 8 lines not only deletes the row for the user but also avoids an error.
        try:
            for item in items_selected:                                            # for each item in the list...
                if item.column() == 0 or item.column() == 1 or item.column() == 2: # if the item is a well, replicate, or sample ID...
                    self.settable.removeRow(item.row())                            # Remove the row
                else:
                    self.settable.setItem(item.row(), item.column(), QTableWidgetItem('X')) # else replace the contents with an 'X'
        except RuntimeError:
            print('Deletion Error') # I'm not exactly sure what the error means but this helps prevent the program from crashing

        ### The process is less complicates with the main table
        items_selected = self.maintable.selectedItems()
        for item in items_selected:
            self.maintable.setItem(item.row(), item.column(), QTableWidgetItem('X'))

    ### Deletes the selected row(s) in the main table
    def delete(self):
        try:
            slctrow = self.maintable.selectedItems() # Get a list of selected cells
            rls = slctrow[0].row()                   # Row left selected = the row of the first selected cell
            for r in slctrow:                        # for each cell...
                self.maintable.removeRow(r.row())    # Remove the row
                self.netaddedrows -= 1               # Mark this row removal

            try:
                ### Normally removing a row removes that rows selection.
                # By reselecting the row, the user doesn't have to keep clicking
                # the same row to delete multiple rows
                self.maintable.item(rls, 0).setSelected(True)
            except:
                x = 'Do Nothing' #  If the rls can't be selected, just leave nothing selected
        except IndexError:
            x = 'Do Nothing'

    ### Inserts a new row above the currently selected row
    def insert(self):
        try:
            slctrow = self.maintable.selectedItems()[0].row() # Get the fist selected row (only first to avoid confusion and unintentialal row insertions)
            self.maintable.insertRow(slctrow)                 # Insert the row

            for c in range(5): # for each column in the new row...
                self.maintable.setItem(slctrow, c, QTableWidgetItem('X')) # Put an 'X' in each cell
                ### ^ this allows the new row to easlily be deleted if desired

            self.maintable.item(slctrow+1,0).setSelected(False) # Deselect the original row
            self.maintable.item(slctrow,0).setSelected(True)    # Select the new row

            self.netaddedrows += 1 # Mark the row addition
        except IndexError:
            x = 'Do Nothing' # If no row is selected, do nothing

    ### Removes the most recently added block in the main table
    def remove(self):
        ### Find the number of rows that have data and store the number in the 'rm' variable
        rm = 0
        while str(self.maintable.item(rm, 0)) != 'None':
            rm += 1

        for r in range(self.savedrownums[-1] + 1): # for the number of rows in the last addition...
            self.maintable.removeRow(rm - r)       # Remove each of the rows starting from the last

        del self.savedrownums[-1] # Erase the last addition from memory

    ### Clears the contents of the main table
    def clearmaintab(self):
        self.maintable.clearContents() # Clear the contents of the main table
        self.savedrownums = []         # Forget all additions
        self.netaddedrows = 0          # Reset net number of added rows

    ### Adds the block of information in the set table to the end of the information in the main table
    def addtotable(self):
        rs = 0 # Variable to find end of each column
        totrowskip = sum(self.savedrownums) + self.netaddedrows # Total rows = sum of additions + net rows added

        while str(self.settable.item(rs, 0)) != 'None': # while the cell is not empty (i.e. for each relevant row)...
            for c in range(5): # for each column...
                try:
                    self.maintable.setItem(rs + totrowskip, c, QTableWidgetItem(self.settable.item(rs, c).text())) # Transfer the information from the cell in the set table to the appropriate cell in the maintable
                except:
                    break # If an item can't be added, break the loop
            rs += 1 # Move to next row

        self.savedrownums.append(rs) # Add this addition to the memory of additions

    ### Updates the sample ID when the enzyme or substrate IDs are changed
    def sampleid_update(self):
        self.sampleid = self.enzid_lineedit.text() + '~' + self.subid_lineedit.text() # Form the sample ID

        ### Add the new ID to each well in the set table
        r = 0
        while str(self.settable.item(r, 0)) != 'None':
            self.settable.setItem(r, 2, QTableWidgetItem(self.sampleid))
            r += 1

    ### Transfers the selected wells in the headers table to the set table
    def addwells(self):
        prewells = self.headers_table.selectedItems() # Get a list of the selected wells (wells are not necessarilly in order)
        wells = []                                    # Initialize final list of wells
        for well in prewells:                         # for each selected well
            wells.append(well.text())                 # Add each well text to the final list
        wells.sort()                                  # Sort the wells in alphabetical order (if the wells are numerical, it will sort them numerically)

        ### Add each well to the set table in order
        for well in wells:
            self.settable.setItem(wells.index(well), 0, QTableWidgetItem(well))

    ### Updates the replicate number when it is changed
    def update_repl(self):
        # For each well in the set table, set or replace the replicate number
        r = 0
        while str(self.settable.item(r, 0)) != 'None':
            self.settable.setItem(r, 1, QTableWidgetItem(str(self.repl_spinbox.value())))
            r += 1

    ### Autosets each row when a well is added
    def tableupdate(self, item):
        self.sampleid = self.enzid_lineedit.text() + '~' + self.subid_lineedit.text() # Set sample ID
        if item.column() == 0: # if the item added was in the 'Wells' column...
            self.settable.setItem(item.row(), 1, QTableWidgetItem(str(self.repl_spinbox.value()))) # Set the replicate number
            self.settable.setItem(item.row(), 2, QTableWidgetItem(str(self.sampleid)))             # Set the sample ID
            if str(self.settable.item(item.row(), 3)) == 'None':            # if there is nothing already in the substrate concentration column...
                self.settable.setItem(item.row(), 3, QTableWidgetItem('X')) # Set an 'X'
            if str(self.settable.item(item.row(), 4)) == 'None':            # if there is nothing already in the enzyme concentration column...
                self.settable.setItem(item.row(), 4, QTableWidgetItem('X')) # Set an 'X'

'Application excecution'
def main():
    if hasattr(sys, 'frozen'):
        gs.appdir = sys.executable
        if platform.system() == 'Windows':
            gs.appdir = gs.appdir[:gs.appdir.rfind("\\") + 1]
        else:
            gs.appdir = gs.appdir[:gs.appdir.rfind/("/") + 1]
    else:
        gs.appdir = os.path.dirname(os.path.abspath(__file__))
        if platform.system() == 'Windows':
            gs.appdir += '\\'
        else:
            gs.appdir += '/'

    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    app = Qt.QApplication(sys.argv)
    app.setOrganizationName("TrinhLab-UTK")
    app.setApplicationName("KmCalc-auto formatter")
    gs.rows_errDialog = rows_errDialog()
    gs.DataConf = DataConf()
    gs.colMap = colMap()
    gs.DataSelection = DataSelection(os.getcwd())
    gs.DataSelection.show()
    print("Current app directory is: " + gs.appdir + "\n")
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
