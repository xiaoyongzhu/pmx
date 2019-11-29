#!/usr/bin/env python

import argparse
import logging
import os
import re
import shutil as sh
import sys
from pmx.scripts.cli import check_unknown_cmd


# ==============================================================================
#                             HELPER CLASSES
# ==============================================================================
class NoMissingModuleFilter(logging.Filter):
    def filter(self, record):
        keywords=["module without the python package"]
        return not any(s in record.getMessage() for s in keywords)

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
    """
    Parse mdp file for expected end time and time between saved frames.

    Returns:
    ----------
    end_time: float
        time of last frame (ps).
    dtframe: float
        time between saved frames (ps).
    """
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
                matchObj = re.match( r'tinit\s*=\s*([-+]?[0-9]*\.?[0-9]+)', line, re.M|re.I)
                tinit = float(matchObj.group(1))
            elif("dt" in line):
                matchObj = re.match( r'dt\s*=\s*([-+]?[0-9]*\.?[0-9]+)', line, re.M|re.I)
                dt = float(matchObj.group(1))
            elif("nstxout" in line):
                #print(fname,":\t",line)
                matchObj = re.match( r'nstxout\s*=\s*(\d+)', line, re.M|re.I)
                nstxout = int(matchObj.group(1))
    end_time=tinit+(nsteps*dt)
    dtframe=nstxout*dt
    return(end_time, dtframe)


# raw_input returns the empty string for "enter"
def confirm_defNO(msg):
    yes = {'y', 'yes', 'ye'}
    no = {'n','no', ''}
    sys.stderr.write(msg + "\t(y/N):\n")
    choice = input().lower()
    if choice in yes:
        return True
    elif choice in no:
        return False
    else:
        sys.stderr.write("Unexpected input string. Assuming no.\n")
        return False


# ==============================================================================
#                      COMMAND LINE OPTIONS AND MAIN
# ==============================================================================
def parse_options(SGE=False):
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
            'and the non-equilibrium method.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--toppath',
                        type=str,
                        dest='toppath',
                        help='Path to itp and structure files describing'
                            ' the protein and the ligand.',
                        default='../data')
    parser.add_argument('--basepath',
                        type=str,
                        dest='basepath',
                        help='Path where all everything will be done.',
                        default=os.getcwd())
    parser.add_argument('--mdppath',
                        type=str,
                        dest='mdppath',
                        help='Path to mdp files for'
                            ' the protein and ligand simulations.',
                        default='../data/mdp')
    parser.add_argument('--bt',
                        dest='bt',
                        type=str,
                        choices=['triclinic', 'cubic',
                                 'dodecahedron', 'octahedron'],
                        help='Box type.',
                        default='dodecahedron')
    parser.add_argument('-d',
                        dest='d',
                        type=float,
                        help='Distance (nm) between the solute and the box.',
                        default=1.5)
    parser.add_argument('-b',
                        dest='b',
                        type=float,
                        help='Time (ps) at which to start sampling frames '
                        'from equilibrium simulations.',
                        default=2256.0)
    # parser.add_argument('--gmx',
    #                     dest='gmx',
    #                     type=str,
    #                     help='Call to gmx',
    #                     default="gmx")


    if(SGE):
        parser.add_argument('--pe',
                            dest='pe',
                            type=str,
                            help='Parellel environment to use for SGE jobs.',
                            default="openmp_fast")
        parser.set_defaults(rem_sched=False)
        parser.add_argument('--rem_sched',  action='store_true',
                            dest='rem_sched',
                            help='If supplied, luigi will use a central scheduling '
                            'server to manage tasks in the workflow. '
                            'hostname and port will be read from luigi\'s '
                            'config file. '
                            'Otherwize, will use a local scheduler on the '
                            'current machine.')
        parser.add_argument('--mdrun',
                            dest='mdrun',
                            type=str,
                            help='Call to mdrun. For best performance on a cluster '
                            'make this value aliased for the optimal version of '
                            'mdrun on each node.',
                            default="gmx mdrun")
    else:
        parser.add_argument('--mdrun',
                            dest='mdrun',
                            type=str,
                            help='Call to mdrun.',
                            default="gmx mdrun")

    parser.add_argument('--mdrun_opts',
                        dest='mdrun_opts',
                        type=str,
                        help='Optional arguments to mdrun. '
                        'Enclose in quotes.',
                        default="")



    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args
