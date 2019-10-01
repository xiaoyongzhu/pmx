#!/usr/bin/env python

import argparse
import glob
import numpy as np
import os
import shutil as sh
import sys
import warnings
from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
from pmx.model import Model
from pmx.scripts.cli import check_unknown_cmd

# Constants
kb = 0.00831447215   # kJ/(K*mol)


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
    

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args



# ==============================================================================
#                             Workflow Class
# ==============================================================================
class Workflow_inProtein:
    def __init__(self, toppath, mdppath, proteins=[], ligands=[],
                 n_repeats=3, n_sampling_sims=1, basepath=os.getcwd(),
                 d=1.5, bt="dodecahedron", salt_conc=0.15,
                 mdrun="gmx mdrun", mdrun_opts=""):
        self.toppath = toppath
        self.mdppath = mdppath
        self.n_repeats = n_repeats
        self.n_sampling_sims = n_sampling_sims
        self.proteins = proteins
        self.ligands = ligands
        self.basepath = basepath
        self.d = d
        self.bt = bt
        self.salt_conc = salt_conc
        self.mdrun = mdrun
        self.mdrun_opts = mdrun_opts
        self.states={"A":"l0", "B":"l1"} #states and suffixes of mdp files
        
        
    def check_sanity(self):
        if(not sh.which('gmx')):
            raise RuntimeError('gmx not found in $PATH!')
        if(not sh.which('perl')):
            raise RuntimeError('perl not found in $PATH!')
        
    def check_inputs(self):
        if(not os.path.isdir(self.toppath)):
            raise RuntimeError('Folder %s provided as toppath'
                               'does not exist'%self.toppath)
            
    def gen_folder_name(self,protein,ligand):
        return(self.basepath+'/prot_'+protein+'/lig_'+ligand)
                               
    def gather_inputs(self, clean=True):
        for p in self.proteins:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                os.makedirs(folder, exist_ok=True)
                os.chdir(folder)
                
                #topology
                sh.copy(self.toppath+"/topol_abs_prot_norestr_amber.top",
                        folder+"/topol.top")
                sh.copy(self.toppath+"/ligand/"+l+"/lig.itp",folder+"/")
                sh.copy(self.toppath+"/proteins/"+p+"/prot.itp",folder+"/")
                
                #initial coordinates where protein and ligand are bound
                sh.copy(self.toppath+"/proteins/"+p+"/prot_lig.pdb",
                        folder+"/init.pdb")
                
                #generate temporary index file
                os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                          "-o index.ndx > setup.log 2>&1")
                
                #generate restraints for equillibration
                #TODO: rewrite this to use the pmx Topology class
                os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                          "-fc 9000 9000 9000 -o prot_posre.itp "
                          "-n index.ndx >> setup.log 2>&1")
                check_file_ready("prot_posre.itp")
                os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                          "-fc 500 500 500 -o prot_posre_soft.itp "
                          "-n index.ndx >> setup.log 2>&1")
                check_file_ready("prot_posre_soft.itp")
                os.system("echo 'MOL\n' | gmx editconf -f init.pdb "
                          "-o lig.pdb -n index.ndx >> setup.log 2>&1")
                check_file_ready("index.ndx")
                os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                          "-fc 9000 9000 9000 "
                          "-o lig_posre.itp >> setup.log  2>&1")
                check_file_ready("lig_posre.itp")
                os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                          "-fc 500 500 500 "
                          "-o lig_posre_soft.itp >> setup.log 2>&1")
                check_file_ready("lig_posre_soft.itp")
                
                #clean overwritten files
                if(clean):
                    cleanList = glob.glob(folder+'/#*')
                    for filePath in cleanList:
                        try:
                            os.unlink(filePath)
                        except:
                            print("Error while deleting file: ", filePath)
                
                #Return to basepath
                os.chdir(self.basepath)
                
    def prep(self, clean=True):
        for p in self.proteins:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                os.chdir(folder)
                
                #generate the box and solvate
                if(not os.path.isfile("tpr.tpr")): #skip if it already exists
                    os.system("gmx editconf -f init.pdb -o box.pdb -bt %s -d %f "\
                              "> prep.log 2>&1"%(self.bt, self.d))
                    check_file_ready("box.pdb")
                    os.system("gmx solvate -scale 1.0 -cp box.pdb -o water.pdb "\
                              "-cs spc216.gro -p topol.top >> prep.log 2>&1")
                    check_file_ready("water.pdb")
                    os.system("gmx grompp -p topol.top -c water.pdb -o tpr.tpr "\
                              "-f %s/protein/init.mdp -v -maxwarn 2 "\
                              ">> prep.log 2>&1"%self.mdppath)
                    check_file_ready("tpr.tpr")
                
                #generate ions for each
                #independent repeat (multiple for confidence estimate)
                for i in range(self.n_repeats):
                    #sampling simulations in each repeat
                    for m in range(self.n_sampling_sims):
                        top_ions="topol_ions%d_%d.top"%(i,m)
                        pdb_ions="ions%d_%d.pdb"%(i,m)
                        if(os.path.isfile(pdb_ions)): #skip if it already exists
                            continue
                        sh.copy("topol.top", top_ions)
                        os.system("echo 'SOL' | gmx genion -s tpr.tpr "
                                  "-p %s -conc %f "
                                  "-neutral -nname Cl -pname Na "
                                  "-o %s >> genion.log 2>&1" %(
                                      top_ions, self.salt_conc, pdb_ions) )
                        check_file_ready(pdb_ions)
                                        
                #clean overwritten files
                if(clean):
                    cleanList = glob.glob(folder+'/#*')
                    for filePath in cleanList:
                        try:
                            os.unlink(filePath)
                        except:
                            print("Error while deleting file: ", filePath)
                
                #Return to basepath
                os.chdir(self.basepath)
                
    def run_stage(self, stage, mdp_template, struct_template,
                  posre=None, clean=True):
        print("Running stage "+stage+":")
        for p in self.proteins:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                
                #independent repeats for error analysis
                for i in range(self.n_repeats):
                    #sampling simulations in each repeat
                    for m in range(self.n_sampling_sims):
                        for s in self.states:
                            sim_folder=folder+"/state%s/repeat%d/%s%d"%(s,i,stage,m)
                            os.makedirs(sim_folder, exist_ok=True)
                            os.chdir(sim_folder)
                            print("\t"+sim_folder)
                            #don't rerun if already ready
                            if(os.path.isfile("confout.gro")):
                                print("\t\tPrevious")
                                continue
                            
                            #"../data/mdp/em_{0}.mdp"
                            mdp=mdp_template.format(self.states[s])#insert the l0/l1 suffix
                            #"prot_{0}/lig_{1}/ions{3}_{4}.pdb"
                            #"prot_{0}/lig_{1}/state{2}/repeat{3}/em{4}/confout.gro"
                            struct=struct_template.format(p,l,s,i,m)#insert the s,i, and m suffixes
                            top_template=folder+"/topol_ions{3}_{4}.top"
                            top=top_template.format(p,l,s,i,m)
                            
                            #restraint
                            r=""
                            if posre:
                                r="-r "+posre.format(p,l,s,i,m)
                            
                            #make tpr
                            os.system("gmx grompp -p %s -c %s %s "
                                      "-o tpr.tpr -f %s -v -maxwarn 3 "
                                      "> prep.log 2>&1"%(
                                          top, struct, r, mdp)
                                      )
                            check_file_ready("tpr.tpr")
                            
                            #run sim
                            os.system(self.mdrun+" -s tpr.tpr "+\
                                      self.mdrun_opts+" > mdrun.log 2>&1")
                            check_file_ready("confout.gro")
                            print("\t\tDone")
                            
                            
                            #clean overwritten files
                            if(clean):
                                cleanList = glob.glob(sim_folder+'/#*')
                                for filePath in cleanList:
                                    try:
                                        os.unlink(filePath)
                                    except:
                                        print("Error while deleting file: ", filePath)
                            
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
    toppath=os.path.abspath(args.toppath)
    mdppath=os.path.abspath(args.mdppath)
    basepath=os.getcwd()
    
    w=Workflow_inProtein(toppath, mdppath, ["BRD1"], ["lig"],
                         basepath=basepath,
                         mdrun="mdrun_threads_AVX2_256",
                         mdrun_opts="-pin on -nsteps 10")
    
    #sanity checks
    w.check_sanity()
    w.check_inputs()
        
    #copy data (*.itp, template topology, ligand and protein structures) to CWD
    w.gather_inputs()
    
    #solvate and generate ions
    w.prep()
    
    #run EM
    w.run_stage("em", mdppath+"/protein/em_posre_{0}.mdp",
                basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                posre=basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb")
        
    #run NVT w hard position restraints to prevent protein deformation
    w.run_stage("nvt_posre", mdppath+"/protein/eq_nvt_posre_{0}.mdp",
                basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/em{4}/confout.gro",
                posre=basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb")
    
    #run NVT w softer position restraints
    w.run_stage("nvt_posre_soft", mdppath+"/protein/eq_nvt_posre_soft_{0}.mdp",
                basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre{4}/confout.gro",
                posre=basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb")
    
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
