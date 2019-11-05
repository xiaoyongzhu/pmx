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

# ==============================================================================
#                             Workflow Class
# ==============================================================================
class Workflow:
    def __init__(self, toppath, mdppath, hosts=[], ligands=[],
                 n_repeats=3, n_sampling_sims=1, basepath=os.getcwd(),
                 d=1.5, bt="dodecahedron", salt_conc=0.15,
                 mdrun="gmx mdrun", mdrun_opts=""):
        self.toppath = toppath
        self.mdppath = mdppath
        self.n_repeats = n_repeats
        self.n_sampling_sims = n_sampling_sims
        self.hosts = hosts
        self.ligands = ligands
        self.basepath = basepath
        self.d = d
        self.bt = bt
        self.salt_conc = salt_conc
        self.mdrun = mdrun
        self.mdrun_opts = mdrun_opts
        self.states=[]
        
        
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
                                               
    def run_callback_on_folders(self, stage, callbackfunc,
                                completition_check=None, **kwargs):
        print("Running stage "+stage+":")
        for p in self.hosts:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                
                if(completition_check and
                   os.path.isfile(folder+"/"+completition_check)):
                        print("\t\tPrevious")
                        continue
                # mydict={'folder':folder,'p':p,'l':l,
                #         'states':self.states,
                #         'n_repeats':self.n_repeats,
                #         'n_sampling_sims':self.n_sampling_sims,
                #         'stage':stage,
                #         'basepath':self.basepath}
                mydict={'folder':folder,'p':p,'l':l,'stage':stage}
                kwargs.update(mydict)
                callbackfunc(**kwargs)
                print("\t\tDone")

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
                         mdrun_opts="-pin on -nsteps 1000")
    
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
                posre=basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                completition_check="confout.gro")
        
    #run NVT w hard position restraints to prevent protein deformation
    w.run_stage("nvt_posre", mdppath+"/protein/eq_nvt_posre_{0}.mdp",
                basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/em{4}/confout.gro",
                posre=basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                completition_check="confout.gro")
    
    #run NVT w softer position restraints
    w.run_stage("nvt_posre_soft", mdppath+"/protein/eq_nvt_posre_soft_{0}.mdp",
                basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre{4}/confout.gro",
                posre=basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                completition_check="confout.gro")
    
    #run NPT to sample starting frames for TI
    w.run_stage("npt", mdppath+"/protein/eq_npt_test_{0}.mdp",
                basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre_soft{4}/confout.gro",
                completition_check="confout.gro")
    
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
