#!/usr/bin/env python

import argparse
import glob
import numpy as np
import os
import shutil as sh
from pmx import ndx, geometry
from pmx.analysis import read_dgdl_files, plot_work_dist, ks_norm_test
from pmx.model import Model
from pmx.scripts.cli import check_unknown_cmd
from pmx.xtc import Trajectory
#from pmx.scripts.workflows.Workflow import Workflow, check_file_ready
from Workflow import Workflow, check_file_ready, copy_if_missing
from find_avg import find_avg_struct
from find_anchors_and_write_ii import find_restraints
from fit_ligs_multiframes_python3 import fit,rotate_velocities_R

# Constants
kb = 0.00831447215   # kJ/(K*mol)


# ==============================================================================
#                            HELPER FUNCTIONS
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
    

    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args



# ==============================================================================
#                             Workflow Class
# ==============================================================================
class Workflow_alligned_inProtein(Workflow):
    def __init__(self, toppath, mdppath, hosts=[], ligands=[],
                 n_repeats=3, n_sampling_sims=1, basepath=os.getcwd(),
                 d=1.5, bt="dodecahedron", salt_conc=0.15,
                 mdrun="gmx mdrun", mdrun_opts=""):
        Workflow.__init__(self, toppath, mdppath, hosts, ligands,
                          n_repeats, n_sampling_sims, basepath,
                          d, bt, salt_conc, mdrun, mdrun_opts) 
        self.states={"A":"l0", "B":"l1"} #states and suffixes of mdp files
        self.TIstates={"A":"l0", "C":"l1"} #states and suffixes of mdp files
        
        
    def check_sanity(self):
        if(not sh.which('gmx')):
            raise RuntimeError('gmx not found in $PATH!')
        if(not sh.which('perl')):
            raise RuntimeError('perl not found in $PATH!')
        
    def check_inputs(self):
        if(not os.path.isdir(self.toppath)):
            raise RuntimeError('Folder %s provided as toppath'
                               'does not exist'%self.toppath)
            
    def gen_folder_name(self,host,ligand):
        return(self.basepath+'/prot_'+host+'/lig_'+ligand)
                               
    def gather_inputs(self, clean=True):
        for p in self.hosts:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                os.makedirs(folder, exist_ok=True)
                os.chdir(folder)
                
                #topology
                copy_if_missing(self.toppath+"/topol_abs_prot_norestr_amber.top",
                        folder+"/topol.top")
                copy_if_missing(self.toppath+"/ligand/"+l+"/lig.itp",folder+"/lig.itp")
                copy_if_missing(self.toppath+"/proteins/"+p+"/prot.itp",folder+"/prot.itp")
                
                #initial coordinates where protein and ligand are bound
                copy_if_missing(self.toppath+"/proteins/"+p+"/prot_lig.pdb",
                        folder+"/init.pdb")
                
                #generate temporary index file
                if(not os.path.isfile("index.ndx")):
                    os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                              "-o index.ndx > setup.log 2>&1")
                    check_file_ready("index.ndx")
                
                #generate restraints for equillibration
                #TODO: rewrite this to use the pmx Topology class
                if(not os.path.isfile("prot_posre.itp")):
                    os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                              "-fc 9000 9000 9000 -o prot_posre.itp "
                              "-n index.ndx >> setup.log 2>&1")
                    check_file_ready("prot_posre.itp")
                if(not os.path.isfile("prot_posre_soft.itp")):
                    os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                              "-fc 500 500 500 -o prot_posre_soft.itp "
                              "-n index.ndx >> setup.log 2>&1")
                    check_file_ready("prot_posre_soft.itp")
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
                  posre=None, clean=True, completition_check=None):
        print("Running stage "+stage+":")
        for p in self.hosts:
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
                            if(completition_check and
                               os.path.isfile(completition_check)):
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
            
                
    def run_TI(self, stage, mdp_template, struct_template, n_morphs,
                  posre=None, clean=True, completition_check=None,
                  runfolder=None):

        if(not runfolder):
            runfolder=stage
            
        print("Running stage "+stage+":")
        for p in self.hosts:
            for l in self.ligands:
                folder = self.gen_folder_name(p,l)
                
                #independent repeats for error analysis
                for i in range(self.n_repeats):
                    #sampling simulations in each repeat
                    for m in range(self.n_sampling_sims):
                        for s in self.TIstates:
                            sim_folder=folder+"/state%s/repeat%d/%s%d"%(s,i,runfolder,m)
                            os.makedirs(sim_folder, exist_ok=True)
                            os.chdir(sim_folder)
                            print("\t"+sim_folder)
                            
                            for o in range(n_morphs):
                                #don't rerun if already ready
                                if(completition_check and
                                   #"prot_{0}/lig_{1}/repeat{3}/GenMorphs{4}/dHdl{5}.xvg"
                                   os.path.isfile(completition_check.format(p,l,s,i,m,o))):
                                    print("%dp\t"%o, end = '', flush=True)
                                    continue
                                
                                #"../data/mdp/em_{0}.mdp"
                                mdp=mdp_template.format(self.TIstates[s])#insert the l0/l1 suffix
                                #"prot_{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/frame{5}.gro"
                                struct=struct_template.format(p,l,s,i,m,o)#insert the s,i, and m suffixes
                                top_template=folder+"/topol_ions{3}_{4}.top"
                                top=top_template.format(p,l,s,i,m)
                                
                                #restraint
                                r=""
                                if posre:
                                    r="-r "+posre.format(p,l,s,i,m)
                                
                                #make tpr
                                tprname = "tpr%d.tpr"%o
                                os.system("gmx grompp -p %s -c %s %s "
                                          "-o %s -f %s -v -maxwarn 3 "
                                          "> prep%d.log 2>&1"%(
                                              top, struct, r, tprname, mdp, o)
                                          )
                                check_file_ready(tprname)
                                
                                #run sim
                                os.system(self.mdrun+(" -deffnm tpr%d "%o)+
                                          ("-dhdl dHdl%d.xvg "%o)+
                                          self.mdrun_opts+
                                          (" > mdrun%d.log 2>&1"%o)
                                          )
                                check_file_ready("dHdl%d.xvg"%o)
                                check_file_ready("tpr%d.gro"%o)
                                print("%dd\t"%o, end = '', flush=True)
                                
                            print("\n\t\tDone")
                            
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
                            
    def run_everything(self):
        #sanity checks
        self.check_sanity()
        self.check_inputs()
            
        #copy data (*.itp, template topology, ligand and protein structures) to CWD
        self.gather_inputs()
        
        #solvate and generate ions
        self.prep()
        
        #run EM
        self.run_stage("em", self.mdppath+"/protein/em_posre_{0}.mdp",
                    self.basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                    posre=self.basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                    completition_check="confout.gro")
            
        #run NVT w hard position restraints to prevent protein deformation
        self.run_stage("nvt_posre", self.mdppath+"/protein/eq_nvt_posre_{0}.mdp",
                    self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/em{4}/confout.gro",
                    posre=self.basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                    completition_check="confout.gro")
        
        #run NVT w softer position restraints
        self.run_stage("nvt_posre_soft", self.mdppath+"/protein/eq_nvt_posre_soft_{0}.mdp",
                    self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre{4}/confout.gro",
                    posre=self.basepath+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb",
                    completition_check="confout.gro")
        
        #run NPT to sample starting frames for TI
        self.run_stage("npt", self.mdppath+"/protein/eq_npt_test_{0}.mdp",
                    self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre_soft{4}/confout.gro",
                    completition_check="confout.gro")
                
        #genergate Boresh-style protein-ligand restraints
        self.run_callback_on_folders("Restraints", self.gen_restr_callback,
                    srctpr =self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/npt{4}/tpr.tpr",
                    srctraj=self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/npt{4}/traj.trr",
                    mdp=self.mdppath+"/protein/init.mdp",
                    completition_check="restraint_coord_distrib.png",
                    b=0)
                        
        #align vaccum ligand onto apo protein structures
        self.run_callback_on_folders("GenMorphs", self.gen_alligned_morphs_callback,
                    srctpr =self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/npt{4}/tpr.tpr",
                    srctraj=self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/npt{4}/traj.trr",
                    b=0, n_morphs=21)
        
        #run TI
        self.run_TI("TI", self.mdppath+"/protein/ti_{0}.mdp",
                    self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/frame{5}.gro",
                    n_morphs=21,
                    completition_check=self.basepath+"/prot_{0}/lig_{1}/state{2}/repeat{3}/GenMorphs{4}/tpr{5}.gro",
                    runfolder="GenMorphs")
        
        #analyse dHdl files
        
        #plot summary
    
 
# ==============================================================================
#                           CALLBACK FUNCTIONS
# ==============================================================================                           
        
    def gen_restr_callback(self, **kwargs):
        """Generates Boresh-style restraints by fitting a
        Gaussians to the restraint coordinate distributions.
        
        Selects anchor atoms to maximize Gaussianity of
        restraint distributions and maximise the
        correlation between bound ligand and apo protein.
        This is usefull when the binding pocket
        changes shape on binding/unbinding.
        
        Please use as a callback function by providing it
        to Workflow.run_callback_on_folders()
    
        Parameters
        ----------
        **kwargs: dict        
            Generated by Workflow.run_callback_on_folders()   
            
            b: float
                start time (in ps) for fit.
            srctraj: str
                Source trajectory path string with format keys for protein,
                ligand, state, repeat, and sim id. Eg.:
                    "/prot_{0}/lig_{1}/state{2}/repeat{3}/npt{4}/traj.trr"
            srctpr: str
                Source tpr path string. Same format as srctraj.
            mdp: str
                mdp file that describes the system. Used to make a tpr
                without solvent/ions.
    
        Returns
        -------
        None.
    
        """
        folder = kwargs.get('folder')
        p = kwargs.get('p')
        l = kwargs.get('l')
        
        b=kwargs.get('b', 0) #begining of trj to use (ps)
        srctraj=kwargs.get('srctraj')
        srctpr=kwargs.get('srctpr')
        mdp=kwargs.get('mdp')
        
    
        print(folder)
        os.chdir(folder)
        
        #create prot+MOL index group
        os.system("echo \"1|13\nq\n\" | "
                  "gmx make_ndx -f ions0_0.pdb "
                  "-o index_prot_mol.ndx > /dev/null 2>&1")
        check_file_ready("index_prot_mol.ndx")
                  
        #make topology
        os.system("sed 's/SOL/;SOL/g' topol.top > topol_prot_mol.top")
                          
        for s in self.states:
            #make tprs
            if(s == "A"):   #align A to initial structure
                ref="box.pdb"
            else:           #align B to average of A
                ref="averageA.gro"
            os.system("gmx grompp -p topol_prot_mol.top -c %s "
                      "-f %s -o tpr%s.tpr "
                      "-maxwarn 2 > grompp%s.log 2>&1"%(
                              ref, mdp, s,s) )
            check_file_ready("tpr%s.tpr"%s)
              
            #collect trjs
            if(not os.path.isfile("all_eq%s_fit.xtc"%s)):
                print("\tCollecting trajectories for state%s"%s)
                
                #independent repeats for error analysis
                for i in range(self.n_repeats):
                    #sampling simulations in each repeat
                    for m in range(self.n_sampling_sims):
                        if(not os.path.isfile("eq%s%d_%d.xtc"%(s,i,m))):
                            tpr=srctpr.format(p,l,s,i,m)
                            trj=srctraj.format(p,l,s,i,m)
                            os.system("echo 4 Protein_MOL | "
                                      "gmx trjconv -s %s -f %s "
                                      "-o eq%s%d_%d.xtc "
                                      "-sep -ur compact -pbc mol -center "
                                      "-boxcenter zero -n index_prot_mol.ndx "
                                      "-b %d > /dev/null 2>&1"%(
                                              tpr,trj, s,i,m, b) )
                            check_file_ready("eq%s%d_%d.xtc"%(s,i,m))
                            
                #concatenate trajectories
                os.system("gmx trjcat -f eq%s*.xtc -o all_eq%s.xtc -sort "
                          "-cat > /dev/null 2>&1"%(s,s) )
                check_file_ready("all_eq%s.xtc"%s)
                          
                #fit to reference structure in tpr files
                os.system("echo 4 0 | gmx trjconv -s tpr%s.tpr -f all_eq%s.xtc "
                          "-o all_eq%s_fit.xtc -fit rot+trans > /dev/null 2>&1"%(
                                  s,s,s) )
                check_file_ready("all_eq%s_fit.xtc"%s)
            
            #dump first frame
            if(not os.path.isfile("dump%s.gro"%s)):
                os.system("echo 0 | gmx trjconv -f all_eq%s_fit.xtc "
                          "-s tpr%s.tpr -o dump%s.gro "
                          "-dump 0 > /dev/null 2>&1"%(
                                  s,s,s) )
                check_file_ready("dump%s.gro"%s)
                
            #find avg structure of A
            if(s=="A"):
                if(not os.path.isfile("averageA.gro")):
                    find_avg_struct("dumpA.gro", "all_eqA_fit.xtc",
                                    "averageA.gro")
                    check_file_ready("averageA.gro")
                    
        #generate the restraints
        if(not os.path.isfile("ii.itp")):
            find_restraints(log=False)
            check_file_ready("ii.itp")
            check_file_ready("out_dg.dat")
        
        #restore base path    
        os.chdir(self.basepath)
                            
        
                            
    def gen_alligned_morphs_callback(self, **kwargs):
        """Generates coordinates for bound restrained but decoupled state (C)
        by structure allignemnt of vacuum ligand into apo protein.
        Initial frames for TI in coupled state (A) are generated by
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
        
        #independent repeats for error analysis
        for i in range(self.n_repeats):
            #sampling simulations in each repeat
            for m in range(self.n_sampling_sims):
                #generate morphs for A state
                s="A"
                frname="frame.gro"
                
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
                
                os.chdir(self.basepath)
                    
                #now make the C state
                s="C"
                sim_folder=folder+"/state%s/repeat%d/%s%d"%(s,i,stage,m)
                os.makedirs(sim_folder, exist_ok=True)
                os.chdir(sim_folder)
                
                ready=False
                for o in range(n_morphs):
                    if(os.path.isfile("frame%d.gro"%o)):
                        ready=(os.path.getctime(trj) <
                               os.path.getctime('frame0.gro'))
                        if(not ready):
                            break
                
                if(not ready): # at least one is missing
                    m_A = Model(folder+"/ions%d_%d.pdb"%(i,m),bPDBTER=True)
                    m_B = Model(folder+"/ions%d_%d.pdb"%(i,m),bPDBTER=True)
                    m_C = Model(folder+"/ions%d_%d.pdb"%(i,m),bPDBTER=True)
                    m_A.a2nm()
                    m_B.a2nm()
                    m_C.a2nm()
                    
                    trj_A = Trajectory(srctraj.format(p,l,"A",i,m))
                    trj_B  = Trajectory(srctraj.format(p,l,"B",i,m))
                    
                    ndx_file = ndx.IndexFile(folder+"/index_prot_mol.ndx", verbose=False)
                    p_ndx = np.asarray(ndx_file["Protein"].ids)-1
                    l_ndx = np.asarray(ndx_file["MOL"].ids)-1
                    
                    #Frames are not acessible individually, just in sequence.
                    #pmx.xtc.Trajectory is based on __iter__, so we need a custom
                    #"for" loop to simultaneously go through both trajectories.
                    #Based on https://www.programiz.com/python-programming/iterator
                    iter_A = iter(trj_A)
                    iter_B = iter(trj_B)
                    fridx=0
                    while True:
                        try:
                            frame_A = next(iter_A)
                            frame_B = next(iter_B)
                        except StopIteration:
                            break
                        
                        if(not os.path.isfile("frame%d.gro"%fridx)):
                            frame_A.update(m_A)
                            frame_B.update(m_B)
                            frame_B.update(m_C)
                            
                            # step1: fit prot from prot+lig onto apo protein
                            (v1,v2,R) = fit( m_B, m_A, p_ndx, p_ndx )
                            # rotate velocities
                            # not needed. We aren't saving m_A
                                            
                            # step2: ligand onto the ligand from prot+lig structure
                            (v1,v2,R) = fit( m_A, m_C, l_ndx, l_ndx )
                            # rotate velocities
                            rotate_velocities_R( m_C, R )
                            
                            #replace coordinates and velocities of ligand in B with rotated ones from C
                            for i in l_ndx:
                                m_B.atoms[i].x = m_C.atoms[i].x
                                m_B.atoms[i].v = m_C.atoms[i].v
                    
                            # output                        
                            m_B.write("frame%d.gro"%fridx)
                            
                        fridx+=1
    
                
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
    
    w=Workflow_alligned_inProtein(toppath, mdppath, ["BRD1"], ["lig"],
                         basepath=basepath,
                         #mdrun="mdrun_threads_AVX2_256",
                         mdrun="gmx mdrun",
                         mdrun_opts="-pin on -nsteps 1000 -ntomp 8")

    w.run_everything()

    print("Complete.\n")

if __name__ == '__main__':
    args = parse_options()
    main(args)
