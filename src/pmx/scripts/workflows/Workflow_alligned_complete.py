#!/usr/bin/env python

import os
#from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
from pmx.scripts.cli import check_unknown_cmd
#from pmx.scripts.workflows.Workflow import Workflow, check_file_ready
from Workflow import Workflow
from Workflow_alligned_in_protein import Workflow_alligned_inProtein, parse_options
from Workflow_in_water import Workflow_inWater



# Constants
kb = 0.00831447215   # kJ/(K*mol)


# ==============================================================================
#                             Workflow Class
# ==============================================================================
class Workflow_alligned_complete(Workflow):
    def __init__(self, toppath, mdppath, hosts=[], ligands=[],
                 n_repeats=3, n_sampling_sims=1, basepath=os.getcwd(),
                 d=1.5, bt="dodecahedron", salt_conc=0.15,
                 mdrun="gmx mdrun", mdrun_opts=""):
        Workflow.__init__(self, toppath, mdppath, hosts, ligands,
                          n_repeats, n_sampling_sims, basepath,
                          d, bt, salt_conc, mdrun, mdrun_opts) 

        self.wp = Workflow_alligned_inProtein(toppath, mdppath, hosts, ligands,
                                   n_repeats, n_sampling_sims, basepath,
                                   d, bt, salt_conc, mdrun, mdrun_opts)
        self.ww = Workflow_inWater(toppath, mdppath, ["water"], ligands,
                                   n_repeats, n_sampling_sims, basepath,
                                   d, bt, salt_conc, mdrun, mdrun_opts)
        
    
                            
    def run_everything(self):
        """Runs the whole workflow.
        
        Parameters
        ----------
        None.
    
        Returns
        -------
        None.
        """
        
        #run the lig+protein in water
        self.wp.run_everything()
        
        #run the lig in water
        self.ww.run_everything()
        
        #analyse dHdl files
        
        #plot summary

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
    basepath=os.path.abspath(args.basepath)
    
    w=Workflow_alligned_complete(toppath, mdppath, ["BRD1"], ["lig"],
                         basepath=basepath,
                         mdrun="mdrun_threads_AVX2_256",
                         #mdrun="gmx mdrun",
                         mdrun_opts="-pin on -nsteps 1000 -ntomp 8")

    w.run_everything()

    print("Complete.\n")

if __name__ == '__main__':
    args = parse_options()
    main(args)
