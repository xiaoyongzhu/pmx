#!/usr/bin/env python

from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
from pmx import __version__
import sys
import os
import numpy as np
import argparse
import warnings
from pmx.scripts.cli import check_unknown_cmd

# Constants
kb = 0.00831447215   # kJ/(K*mol)


# ==============================================================================
#                               FUNCTIONS
# ==============================================================================



# ==============================================================================
#                      COMMAND LINE OPTIONS AND MAIN
# ==============================================================================
def parse_options():
    """Parse cmd-line options.

    Returns:
    ----------
    args: argparse.Namespace
        The processed command line arguments.
    """

    parser = argparse.ArgumentParser(description='Runs the whole workflow '\
            'for calculation of free energy of ligand binding '\
            'using a single instance of the ligand per box, '\
            'optimized Boresh-style restraints, '\
            'and non-equilibrium .')

    parser.add_argument('-toppath',
                        type=str,
                        dest='toppath',
                        help='Path to itp and structure files describing the protein and the ligand.',
                        default=['../topologies'])
    parser.add_argument('-t',
                        metavar='temperature',
                        dest='temperature',
                        type=float,
                        help='Temperature in Kelvin. Default is 298.15.',
                        default=298.15)
    

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args


# ==============================================================================
#                               FUNCTIONS
# ==============================================================================
def main(args):
    """Run the main script.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments
    """

   #sanity checks
   
   #copy data (*.itp, template topology, ligand and protein structures) to CWD
   
   #solvate and generate ions
   
   #run EM
   
   #run NVT w hard position restraints to prevent protein deformation
   
   #run NVT w softer position restraints
   
   #run NPT to sample starting frames for TI
   
   #genergate Boresh-style protein-ligand restraints
   
   #align vaccum ligand onto apo protein structures
   
   #run TI
   
   #analyse dHdl files
   
   #plot summary



def entry_point():
    args = parse_options()
    main(args)

if __name__ == '__main__':
    entry_point()
