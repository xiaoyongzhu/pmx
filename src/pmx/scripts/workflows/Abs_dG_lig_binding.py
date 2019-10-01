#!/usr/bin/env python

from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
from pmx import __version__
import sys
import os
import numpy as np
import argparse
import warnings
import shutil as sh
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
                        help='Path to itp and structure files describing'
                            ' the protein and the ligand.',
                        default='../data')
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
#                             Workflow Class
# ==============================================================================
class Workflow_inProtein:
    def __init__(self, args, proteins=[], ligands=[], basepath=os.getcwd()):
        self.args = args
        self.proteins = proteins
        self.ligands = ligands
        self.basepath = basepath
        
    def check_sanity(self):
        if(not sh.which('gmx')):
            raise RuntimeError('gmx not found in $PATH!')
        if(not sh.which('perl')):
            raise RuntimeError('perl not found in $PATH!')
        
    def check_inputs(self):
        if(not os.path.isdir(self.args.toppath)):
            raise RuntimeError('Folder \"%s\" provided as toppath does not exist'\
                               %self.args.toppath)
            
    def gen_folder_name(self,protein,ligand):
        return(self.basepath+'/prot_'+protein+'/lig_'+ligand)
                               
    def gather_inputs(self):
        for p in self.proteins:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                os.makedirs(folder, exist_ok=True)
                
                #topology
                sh.copy(self.args.toppath+"/topol_abs_prot_norestr_amber.top",\
                        folder+"/topol.top")
                sh.copy(self.args.toppath+"/ligand/"+l+"/lig.itp",folder+"/")
                sh.copy(self.args.toppath+"/proteins/"+p+"/prot.itp",folder+"/")
                
                #initial coordinates where protein and ligand are bound
                sh.copy(self.args.toppath+"/proteins/"+p+"/prot_lig.pdb",\
                        folder+"/init.pdb")
                
                os.chdir(folder)
                
                #generate temporary index file
                os.system("echo 'q\\n' | "+\
                          "gmx make_ndx -f init.pdb -o index.ndx > setup.log")
                
                #generate restraints for equillibration
                os.system("echo 'Protein\\n' | "+\
                          "gmx genrestr -f init.pdb "+\
                          "-fc 9000 9000 9000 "+\
                          "-o prot_posre.itp -n index.ndx >> setup.log")
                os.system("echo 'MOL\\n' | "+\
                          "gmx genrestr -f init.pdb "+\
                          "-fc 9000 9000 9000 "+\
                          "-o lig_posre.itp -n index.ndx >> setup.log")
                os.system("echo 'Protein\\n' | "+\
                          "gmx genrestr -f init.pdb "+\
                          "-fc 500 500 500 "+\
                          "-o prot_posre_soft.itp -n index.ndx >> setup.log")
                os.system("echo 'MOL\\n' | "+\
                          "gmx genrestr -f init.pdb "+\
                          "-fc 500 500 500 "+\
                          "-o lig_posre_soft.itp -n index.ndx >> setup.log")
                
                #Return to basepath
                os.chdir(self.basepath)

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

    w=Workflow_inProtein(args, ["BRD1"], ["lig"], os.getcwd())
    
    #sanity checks
    w.check_sanity()
    w.check_inputs()
        
    #copy data (*.itp, template topology, ligand and protein structures) to CWD
    w.gather_inputs()
    
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
