#!/usr/bin/env python

import glob
import os
import shutil as sh
from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
#from pmx.scripts.workflows.Workflow import Workflow, check_file_ready
from Workflow import Workflow, check_file_ready, copy_if_missing
from Workflow_alligned_in_protein import Workflow_alligned_inProtein, parse_options

# Constants
kb = 0.00831447215   # kJ/(K*mol)


# ==============================================================================
#                             Workflow Class
# ==============================================================================
class Workflow_inWater(Workflow_alligned_inProtein):
    def __init__(self, toppath, mdppath, hosts=["water"], ligands=[],
                 n_repeats=3, n_sampling_sims=1, basepath=os.getcwd(),
                 d=1.5, bt="dodecahedron", salt_conc=0.15,
                 mdrun="gmx mdrun", mdrun_opts=""):
        Workflow.__init__(self, toppath, mdppath, hosts, ligands,
                          n_repeats, n_sampling_sims, basepath,
                          d, bt, salt_conc, mdrun, mdrun_opts) 
        self.states={"A":"l0", "B":"l1"} #states and suffixes of mdp files
        self.TIstates=self.states #states and suffixes of mdp files
        
    def gen_folder_name(self,host,ligand):
        return(self.basepath+'/'+host+'/lig_'+ligand)
                               
    def gather_inputs(self, clean=True):
        for p in self.hosts:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                os.makedirs(folder, exist_ok=True)
                os.chdir(folder)
                
                #topology
                copy_if_missing(self.toppath+"/topol_abs_water_amber.top",
                        folder+"/topol.top")
                copy_if_missing(self.toppath+"/ligand/"+l+"/lig.itp",folder+"/lig.itp")
                
                #initial coordinates where protein and ligand are bound
                copy_if_missing(self.toppath+"/ligand/"+l+"/ligand.pdb",
                        folder+"/init.pdb")
                
                #generate temporary index file
                if(not os.path.isfile("index.ndx")):
                    os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                              "-o index.ndx > setup.log 2>&1")
                    check_file_ready("index.ndx")
                
                #generate restraints for equillibration
                #TODO: rewrite this to use the pmx Topology class
                if(not os.path.isfile("lig.pdb")):
                    os.system("echo 'MOL\n' | gmx editconf -f init.pdb "
                              "-o lig.pdb -n index.ndx >> setup.log 2>&1")
                    check_file_ready("lig.pdb")
                if(not os.path.isfile("lig_posre.itp")):
                    os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                              "-fc 9000 9000 9000 "
                              "-o lig_posre.itp >> setup.log  2>&1")
                    check_file_ready("lig_posre.itp")
                if(not os.path.isfile("lig_posre_soft.itp")):
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
        for p in self.hosts:
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
                              "-f %s/water/init.mdp -v -maxwarn 2 "\
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
                
    def run_everything(self):
        """Runs the whole workflow.
        
        Parameters
        ----------
        None.
    
        Returns
        -------
        None.
        """
        
        #sanity checks
        self.check_sanity()
        self.check_inputs()
            
        #copy data (*.itp, template topology, ligand and protein structures) to CWD
        self.gather_inputs()
        
        #solvate and generate ions
        self.prep()
        
        #run EM
        self.run_stage("em", self.mdppath+"/water/em_{0}.mdp",
                    self.basepath+"/{0}/lig_{1}/ions{3}_{4}.pdb",
                    posre=self.basepath+"/{0}/lig_{1}/ions{3}_{4}.pdb",
                    completition_check="confout.gro")
            
        #run NVT w hard position restraints to prevent protein deformation
        self.run_stage("nvt", self.mdppath+"/water/eq_nvt_{0}.mdp",
                    self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/em{4}/confout.gro",
                    posre=self.basepath+"/{0}/lig_{1}/ions{3}_{4}.pdb",
                    completition_check="confout.gro")
            
        #run NPT to sample starting frames for TI
        self.run_stage("npt", self.mdppath+"/water/eq_npt_test_{0}.mdp",
                    self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/nvt{4}/confout.gro",
                    completition_check="confout.gro")
        
        #gen morphs
        self.run_callback_on_folders("GenMorphs", self.gen_morphs_callback,
                    srctpr =self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/npt{4}/tpr.tpr",
                    srctraj=self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/npt{4}/traj.trr",
                    b=0, n_morphs=21)
        
        #run TI
        self.run_TI("TI", self.mdppath+"/water/ti_{0}.mdp",
                    self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/frame{5}.gro",
                    n_morphs=21,
                    completition_check=self.basepath+"/{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/tpr{5}.gro",
                    runfolder="GenMorphs")
        
        
# ==============================================================================
#                           CALLBACK FUNCTIONS
# ==============================================================================   
    def gen_morphs_callback(self, **kwargs):
        """Generates coordinates initial frames for TI by
        extracting them from NPT simulations.
        
        Assumes existing equilibrium simulations for states
        A (coupled, unrestrained) and B (decoupled, unrestrained)
        
        Please use as a callback function by providing it
        to Workflow.run_callback_on_folders()
    
        Parameters
        ----------
        b: float
            start time (in ps) for frame extraction.
        srctraj: str
            Source trajectory path string with format keys for protein,
            ligand, state, repeat, and sim id. Eg.:
                "/prot_{0}/lig_{1}/state{2}/repeat{3}/npt{4}/traj.trr"
        srctpr: str
            Source tpr path string. Same format as srctraj.
        **kwargs: dict        
            Generated by Workflow.run_callback_on_folders()
    
        Returns
        -------
        None.
    
        """
        folder = kwargs.get('folder')
        p = kwargs.get('p')
        l = kwargs.get('l')
        stage=kwargs.get('stage','morphs')
        
        b=kwargs.get('b', 0) #begining of trj to use (ps)
        srctraj=kwargs.get('srctraj')
        srctpr=kwargs.get('srctpr')
        n_morphs=kwargs.get('n_morphs')
        
        print("\t"+folder)
        
        frname="frame.gro"
        
        #independent repeats for error analysis
        for i in range(self.n_repeats):
            #sampling simulations in each repeat
            for m in range(self.n_sampling_sims):
                #states
                for s in self.states:
                    sim_folder=folder+"/state%s/repeat%d/%s%d"%(s,i,stage,m)
                    os.makedirs(sim_folder, exist_ok=True)
                    os.chdir(sim_folder)
                    
                    tpr=srctpr.format(p,l,s,i,m)
                    trj=srctraj.format(p,l,s,i,m)
                    
                    #avoid reruning trjconv if inputs aren't new
                    ready=False
                    for o in range(n_morphs):
                        if(os.path.isfile("frame%d.gro"%o)):
                            ready=(os.path.getctime(trj) <
                                   os.path.getctime('frame0.gro'))
                            if(not ready):
                                break
                    if(not ready):
                        #this is slow
                        os.system("echo 0 | gmx trjconv -s %s "
                                  "-f %s -o %s "
                                  "-b %f -sep -ur compact -pbc mol "
                                  "> /dev/null 2>&1"%(tpr,trj,frname,b) )
                    
                    #restore base path    
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
    basepath=os.path.abspath(args.basepath)
    
    w=Workflow_inWater(toppath, mdppath, ["water"], ["lig"],
                         basepath=basepath,
                         #mdrun="mdrun_threads_AVX2_256",
                         mdrun="gmx mdrun",
                         mdrun_opts="-pin on -nsteps 1000 -ntomp 8")
    
    w.run_everything()

    print("Complete.\n")


if __name__ == '__main__':
    args = parse_options()
    main(args)
