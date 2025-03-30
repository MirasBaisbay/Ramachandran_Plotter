""" 
	====================================================================================
	If running script from the command line, functions here are called to parsed user's 
	arguments into the main() function in RamachandranPlotter.py.
	
	Version 2.0.1:
	 - Relies on the easily accessible Biopython package, rather than Phenix as in 
	   versions <2.0
	 - User arguments can be now easily parsed in from the command line (as facilitated 
	   by functions here)
	 - If required, the script could be implemented into existing protein analysis 
	   pipelines by importing this function ( main() ).

	Author information:
	 - Joseph I. J. Ellaway
	 - josephellaway@gmail.com
	 - https://github.com/Joseph-Ellaway
	====================================================================================
"""

import argparse


def CollctUserArgs():
    parser = argparse.ArgumentParser(
        description="Generate averaged Ramachandran plot from PDB files."
    )
    # Modify the pdb argument to accept one or more values
    parser.add_argument(
        '--pdb', '-p', 
        nargs='+',  # Accepts one or more arguments
        required=True,
        help="List of PDB files (separated by spaces)"
    )
    # Other arguments; for example:
    parser.add_argument(
        '--models', '-m',
        type=int,
        default=1,
        help="Number of models to iterate through"
    )
    parser.add_argument(
        '--chains', '-c',
        default='A',
        help="Chain identifier"
    )
    parser.add_argument(
        '--out_dir', '-d',
        default='.',
        help="Output directory"
    )
    parser.add_argument(
        '--plot_type', '-t',
        type=int,
        default=0,
        help="Plot type (integer corresponding to type, e.g. 0 for 'All')"
    )
    parser.add_argument(
        '--file_type', '-f',
        default='png',
        help="Output file type (png or pdf)"
    )
    parser.add_argument(
        '--save', '-s',
        action='store_true',
        help="Flag to save the user data as CSV"
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Enable verbose output"
    )
    
    return parser.parse_args()



def VerboseStatement(verb_boolean, statement):
	"""
	=============================================================================
	If first argument is true, function prints a given statement to command line.
	=============================================================================
	"""

	if verb_boolean: 		# Parsed from argparser
		print(statement)
		# print("\n")
	else:
		pass