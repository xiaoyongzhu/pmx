#!/usr/bin/env python

import os
import shutil as sh


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

