#!/usr/bin/env python

import os
import sys
#from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
import pmx.scripts.analyze_dhdl
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
        
    def run_analysis(self, n_morphs):
        print("Running stage analysis:")
        
        sim_grps=[[self.wp, self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/dHdl{5}.xvg"],
                  [self.ww, self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/dHdl{5}.xvg"]
                  ]
        
        #ligand+protein aq.
        for grp in sim_grps:
            for p in grp[0].hosts:
                for l in self.ligands:
                    folder = grp[0].gen_folder_name(p,l)
                    print("\t"+folder)
                    
                    #independent repeats for error analysis
                    for i in range(self.n_repeats):
                        ana_folder=folder+"/analysis/repeat%d"%i
                        os.makedirs(ana_folder, exist_ok=True)
                        os.chdir(ana_folder)
                        
                        print("%d"%i, end = '', flush=True)
                        
                        if(not os.path.isfile("wplot.png")): #skip if result already exists
                            
                            dHdlA=[]
                            dHdlB=[]
                            dHdlfmt=grp[1]
                            
                            #sampling simulations in each repeat
                            for m in range(self.n_sampling_sims):
                                keys=list(grp[0].TIstates.keys()) #A,B or A,C depending on grp
                                dHdlA.extend([dHdlfmt.format(p,l,keys[0],i,m,o) for o in range(n_morphs)])
                                dHdlB.extend([dHdlfmt.format(p,l,keys[1],i,m,o) for o in range(n_morphs)])
                            
                            orig_argv = sys.argv
                            orig_stdout=sys.stdout
                            orig_stderr=sys.stderr
                            sys.argv = [['analyze_dhdl.py'],\
                                        ['-fA'], dHdlA,
                                        ['-fB'], dHdlB,
                                        ['--nbins', "10"]]
                            sys.argv = [item for sublist in sys.argv for item in sublist]
                            
                            with open("analysis.log", "w") as f:
                                sys.stdout = f
                                sys.stderr = f
                                pmx.scripts.analyze_dhdl.entry_point()
                                
                            sys.argv = orig_argv #reset argv
                            sys.stdout=orig_stdout #reset output
                            sys.stderr=orig_stderr
                            print("d\t", end = '', flush=True)
                            
                        else: # previous
                            print("p\t", end = '', flush=True)
                            
                        os.chdir(self.basepath)#reset cwd                    
                    print() #new line at end of folder
                    
        
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
        self.run_analysis(21)
        
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
