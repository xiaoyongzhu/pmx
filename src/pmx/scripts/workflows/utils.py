#!/usr/bin/env python

import argparse
import os
import re
import shutil as sh
from pmx.scripts.cli import check_unknown_cmd

# ==============================================================================
#                            HELPER FUNCTIONS
# ==============================================================================
def check_file_ready(fname):
    """Checks if a file was sucessfully created
    and gives an informative error if not.

    Parameters
    ----------
    fname : str
        Path to the file beeing checked

    Raises:
    ----------
    OSError:
        If fname is not found.
    """
    if(not os.path.isfile(fname)):
        raise OSError("Failed creating "+fname)

def copy_if_missing(src, trg):
    """Checks if trg file exists and copies it from src if not.

    Parameters
    ----------
    src : str
        Path to source file
    trg : str
        Path to target file

    Raises:
    ----------
    OSError:
        If trg still doesn't exist failed.
    """
    if(not os.path.isfile(trg)):
        sh.copy(src,trg)
        check_file_ready(trg)

def read_from_mdp(fname):
    nsteps = 0
    tinit = 0.0
    dt = 0.0
    nstxout = 0
    with open(fname,"r") as f:
        lineList = f.readlines()
        for line in lineList:
            if("nsteps" in line):
                matchObj = re.match( r'nsteps\s*=\s*(\d+)', line, re.M|re.I)
                nsteps = int(matchObj.group(1))
            elif("tinit" in line):
                matchObj = re.match( r'tinit\s*=\s*(\d+)', line, re.M|re.I)
                tinit = float(matchObj.group(1))
            elif("dt" in line):
                matchObj = re.match( r'dt\s*=\s*(\d+)', line, re.M|re.I)
                dt = float(matchObj.group(1))
            elif("nstxout" in line):
                #print(fname,":\t",line)
                matchObj = re.match( r'nstxout\s*=\s*(\d+)', line, re.M|re.I)
                nstxout = int(matchObj.group(1))
    end_time=tinit+(nsteps*dt)
    return(end_time, nstxout)


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
                        help='Path to itp and structure files describing'
                            ' the protein and the ligand.',
                        default='../data')
    parser.add_argument('-basepath',
                        type=str,
                        dest='basepath',
                        help='Path where all everything will be done.',
                        default=os.getcwd())
    parser.add_argument('-mdppath',
                        type=str,
                        dest='mdppath',
                        help='Path to mdp files for'
                            ' the protein and ligand simulations.',
                        default='../data/mdp')
    parser.add_argument('-t',
                        metavar='temperature',
                        dest='temperature',
                        type=float,
                        help='Temperature in Kelvin.',
                        default=298.15)
    parser.add_argument('-bt',
                        dest='bt',
                        type=str,
                        choices=['triclinic', 'cubic',
                                 'dodecahedron', 'octahedron'],
                        help='Box type',
                        default='dodecahedron')
    parser.add_argument('-d',
                        dest='d',
                        type=float,
                        help='Distance (nm) between the solute and the box',
                        default=1.5)
    parser.add_argument('-b',
                        dest='b',
                        type=float,
                        help='Time (ps) at which to start sampling frames '
                        'from equilibrium simulations',
                        default=2256.0)


    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args
