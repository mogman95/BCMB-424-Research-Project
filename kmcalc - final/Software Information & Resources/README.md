# KmCalc
_KmCalc_ consists of two programs accessed by and connected through a graphical user interface (GUI). One program is used to clean raw microplate-reader enzyme assays data and format it into an output file that can be analyzed by the second program to determine kinetic parameters.

## Installation
1. To acquire the program files:
	
	* Clone the repository with the following command (from [Git Bash in Windows](https://git-scm.com/downloads)
or the terminal in Linux):

	    ~~~
	    git clone https://github.com/mogman95/BCMB-424-Research-Project.git
	    ~~~

## Usage
1. Make sure that you have python3 and all package requirements (found in `requirements.txt`) downloaded before using this program.
2. When using the KmCalc analysis program, make sure your input file is formatted in the way shown in the example files (a template file is provided). The file format must be excel workbook (`.xlsx`).
The only headers that can be modified are the ones containing samples, which in the examples are labeled A1,A2, etc... Conversion factors and units must be consistent with the example (refer to the Section "How do units in _KmCalc_ work?" below for more details). To generate the report file you can use the GUI.

## GUI
1. You can execute the `main.py` file directly to run the program. The home window with the KmCalc logo should pop up. 


## Examples

These can be used to demonstrate the use of the program.

- `Example Data.xlsx`: Demonstration of analysis on real data.
- `Synthetic Data.xlsx`: Used to demonstrate the accuracy of the calculations. Parameters should be K_M = 130 uM | V_max = 0.0137 uM/s | k_cat = 13.7 s^-1
- `Dirty Data.xlsx`: To be used in the autoformatter to demonstrate its utility in data cleaning and formatting.

## How do units in _KmCalc_ work?

KmCalc has a flexible units handling system that will automatically read the units provided in the input file and will convert them to the disired output units. The following units can be handled:
### Time
- Seconds (s)
- Minutes (m)
- Hours (h)
### Concentration
- Molar (M)
- Millimolar (mM)
- Micromolar (uM)
- Nanomolar (nM)
- Picomolar (pM)
### Rate
- Any combination of the above concentration per time units.

## What if the resulting parameters are wrong?

In some rare instances we have observed numerical accuracy problems, this was solved by changing the units of the extinction coefficient ("absorbance2concentration" which is the inverse of the ext. coeff.) and enzyme concentrations into units with higher order of magnitude.
In addition, units that are inappropriately small or large may cause problems such as parameters reaching upper bounds and being equal to zero.

## How does _KmCalc_ work?
### Analysis algorithm:

1. Reads the input file and extracts the data, column map, and conversions. 
2. Identifies the units of the input data.
3. Finds the rates of reaction from the raw data using scipy linear regression.
4. Converts all output data to disired units.
5. Subtracts blank well rate from other rates.
6. Finds kinetic parameters by using scipy non-linear regression to fit data to the Michaelis–Menten equation.
7. Organizes results.
8. Writes and outputs results (if disired).
9. Visualizes data in a plot including averaged data points with standard deviation error bars and fitted curve.

### Auto-formatter workflow:

1. Select raw data file.
2. Select the header row and last row of disired data.
3. Select columns that the user does not want to include in the cleaned data
- Note: blank columns are automatically deleted.
4. Confirm that the selected cleaned data is correct.
5. Create column map:
    1. Select wells that correspond to a certain sample ID (enzyme-substrate combination) and replicate and add them to the staged data editor block by pressing the 'Add selection' button.
    2. Add the enzyme and substrate concentrations corresponding to the added wells.
    3. Once done editing the staged data add it to the finalized data block.
    4. Confirm that the final column map is correct, then press the 'Save Formatted Data' button to output an `.xlsx` file of the cleaned and formatted data ready to be analyzed.

## Help and Credits
KmCalc was originally developed by Sergio Garcia with development and improvement by David Dooley and Sam Cohen all working with Trinh Lab under Dr. Cong Trinh at the University of Tennessee, Knoxville (UTK).

This specific version of KmCalc is being submitted as the final project for the class BCMB 424 at UTK. Previous versions of KmCalc have been provided in the GitHub repository. Any changes between this version ('kmcalc - final') and the 'kmcalc - pre-fall 2022 semester' version were done by Sam Cohen with advising and mentoring from David Dooley.

If you have any questions, comments, concerns, and/or suggestions about KmCalc please contact Sam Cohen (relevant contact information provided below).

Sam Cohen: scohen13@vols.utk.edu

David Dooley: ddooley2@vols.utk.edu