"""
Command line interface
"""

import argparse
import sys
from kmcalc.core import find_parameters


# Command line interface
def cli(args_=None):
    if args_ is None:
        args_ = sys.argv[1:]

    # Parse input
    parser = argparse.ArgumentParser(description='KmCalc')
    parser.add_argument('input_path', help='Path to .xlsx input file')
    parser.add_argument('-o','--output_path', help='Output file path. Default is input_path-output.xlsx', default=None)
    parser.add_argument('-p','--plot', choices=['no','yes'],
                        help='Generate plots of rate vs concentration. They will be stored within the directiory \"plots\" in your output directory', default='yes')
    parser.add_argument('-v','--verbose_level', choices=['no','yes'],
                        help='Level of output as the program executes.', default='yes')
    parser.add_argument('-t','--linear_region_tol',
                        help='Tolerance for absolute values of the second derivative considered to be 0 used to identify the linear region.', default=1.5)

    args = parser.parse_args(args_)
    find_parameters(**vars(args))


if __name__ == "__main__":
    cli()
