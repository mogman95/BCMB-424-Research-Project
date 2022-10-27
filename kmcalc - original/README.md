# KmCalc
_KmCalc_ determines kinetic parameters from raw microplate-reader enzyme assays data through a graphical or command line interface.


If you use _KmCalc_ for your research please cite:

~~~
CITATION_PLACEHOLDER
~~~

## Installation
1. There are two alternative options to acquire the program files:

	* [Download the zip file](https://github.com/TrinhLab/kmcalc/archive/master.zip) and extract it.
	
	* Clone the repository with the following command (from [Git Bash in Windows](https://git-scm.com/downloads)
or the terminal in Linux):

	    ~~~
	    git clone https://github.com/trinhlab/kmcalc.git
	    ~~~

2. Install the python package<sup>[1](#myfootnote1)</sup>:

    ~~~
    pip install -e kmcalc
    ~~~

## Usage
1. Make sure you input file is formatted as [example 1](https://github.com/TrinhLab/kmcalc/blob/master/examples/1/ex1.xlsx) if using `.xlsx` format or [example 2](https://github.com/TrinhLab/kmcalc/blob/master/examples/2/input) for `.csv`.
The only headers that can be modified are the ones containing samples, which in the example are labeled A1,A2, etc... Conversion factors and units must be consistent with the example (refer to the Section "How do units in _KmCalc_ work?" below for more details). To generate the report file you can either use the GUI or CLI:

### GUI
2. Run the command <sup>[1](#myfootnote1)</sup> `kmcalcgui`. If this does not work you can execute the `gui.py` file directly: 

    ~~~
    python kmcalc\kmcalc\gui.py
    ~~~




### CLI
 2. Run the following command <sup>[1](#myfootnote1)</sup>:

    ~~~
    kmcalc <input_file_path>
    ~~~
    For instance, to run the example from the same directory _KmCalc_ was installed, you would execute:

    ~~~
    kmcalc kmcalc/example/ex1.xlsx
    ~~~

To explore additional options execute `kmcalc --help`.

### Examples

These are contained under [kmcalc/examples](https://github.com/TrinhLab/kmcalc/blob/master/examples/):

- Example 1: Simple general usage.
- Example 2: Like example 1 but with `.csv` input format. Note that .csv input may have lower accuracy and lead to a small increase in error bounds of the final parameters.
- Example 3: Illustrates user specified linear region.

## How do units in _KmCalc_ work?

The input and output files include units of conversion factors, enzyme concentrations, and calculated parameters for the sake of clarity. However, the program does not parse this units or performs any sort of unit conversion. Thus if you wish to use different units all you need to do is change the appropriate numeric values. Be aware that the resulting output will still be labeled as having the "default units". In other words, you can consider that all unit specifications in input and output files are for user guidance and that you are free to ignore them and use values in whatever units you wish, so long as they are consistent.

## What if the resulting parameters are wrong?

In some rare instances we have observed numerical accuracy problems, this was solved by changing the units of the extinction coefficient ("absorbance2concentration" which is the inverse of the ext. coeff.) and enzyme concentrations into units with higher order of magnitude.

## How does _KmCalc_ work?
The algorithm operates as follows:

1. Identify the _linear region_ in each data point by examining when the second derivative is close to 0. The second derivative is estimated with a convolution filter which is resistant to noise.
2. Calculate measured product synthesis rate corresponding to V by linear regression of the points in the _linear region_.
3. Subtract blank (substrate concentration is 0) of each replicate to all points of that replicate.
3. Calculate Km and Vmax using robust non-linear regression which is regarded as the most accurate approach [ref](https://www.hindawi.com/journals/jchem/2017/6560983/).
All replicates are included in this calculation which reports parameter variance determined by reduced chi square test [ref](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).
4. Report average and standard deviation for each replicate together with diagnostic information and plots.

## Development
(These topics are not necessary to use _KmCalc_)

### Setup a virtual environment
To setup a virtual environment for _KmCalc_:

~~~
virtualenv -p /usr/bin/python3 /path/to/your/envs/kmcalc
~~~

Gooey depends on wxPython which requires super user privilege to be
installed and thus will not be installed in your environment following the
conventional `pip install -r requirements`
Instead you must install wxpython in your system (for Arch Linux `sudo pacman -S python-wxpython`), then proceed to `pip install -r requirements`, you may
need to copy or simlink the global installation of wxPython to your
environment, but this was not needed in Arch.

### Generate binary executable
Within the virtual environment (this is very important to minimize size) run
 bundle.sh. Prior to this, you may need to copy your global distutil
 installation into the environment, e.g. `cp -r /lib/python-X/distutils /env/kmcalc/lib/python-X/`. This is currently not practical due to the large size of the resulting executable.

## Footnotes
<a name="myfootnote1">1</a>: For Windows is recommended to run commands using the [Anaconda Prompt](https://www.anaconda.com/distribution/), or a Linux emulator such as Cygwin or Windows Subsystem for Linux. More generally look on your search engine of choice the steps needed to run python programs in your operating system of choice.

