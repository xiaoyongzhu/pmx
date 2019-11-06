#!/usr/bin/env python

import glob
import luigi
import os
import sh
from luigi.contrib.sge import LocalSGEJobTask
from workflow import check_file_ready
from Workflow_alligned_in_protein import Workflow_alligned_inProtein
from SGE_tasks.Sim import SGE_Sim

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_PL_folder(LocalSGEJobTask): # will execute on the login node
    base_path = luigi.Parameter(description='Path to the root of dir of the study')
    folder_path = luigi.Parameter(significant=False,
                      description='Path to the protein+ligand folder to set up')
    top_path = luigi.Parameter(significant=False,
                      description='Path to initial topology/structure files')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    
    
    def work(self):
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)
        
        #topology
        sh.copy(self.top_path+"/topol_abs_prot_norestr_amber.top",
                self.folder_path+"/topol.top")
        sh.copy(self.top_path+"/ligand/"+self.l+"/lig.itp",self.folder_path+"/lig.itp")
        sh.copy(self.top_path+"/proteins/"+self.p+"/prot.itp",self.folder_path+"/prot.itp")
                
        #initial coordinates where protein and ligand are bound
        sh.copy(self.top_path+"/proteins/"+self.p+"/prot_lig.pdb",
                self.folder_path+"/init.pdb")
                                
        #generate temporary index file
        os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                  "-o index.ndx > setup.log 2>&1")
        check_file_ready("index.ndx")
        
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
        check_file_ready("lig.pdb")
        os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                  "-fc 9000 9000 9000 "
                  "-o lig_posre.itp >> setup.log  2>&1")
        check_file_ready("lig_posre.itp")
        os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                  "-fc 500 500 500 "
                  "-o lig_posre_soft.itp >> setup.log 2>&1")
        check_file_ready("lig_posre_soft.itp")
        
        #clean overwritten files
        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                raise OSError("Error while deleting file: "+filePath)
        
        #Return to basepath
        os.chdir(self.base_path)
        
    def output(self):
        files=["topol.top", "lig.itp", "prot.itp", "init.pdb", "index.ndx",
               "prot_posre.itp", "prot_posre_soft.itp",
               "lig_posre.itp", "lig_posre_soft.itp"]
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]
    
    
class Prep_PL_folder(LocalSGEJobTask): # will execute on the login node
    base_path = luigi.Parameter(description='Path to the root of dir of the study')
    folder_path = luigi.Parameter(significant=False,
                         description='Path to the protein+ligand folder to set up')
    top_path = luigi.Parameter(significant=False,
                         description='Path to initial topology/structure files')
    mdp_path = luigi.Parameter(significant=False,
                         description='Path to folder with mdp files')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    
    d = luigi.FloatParameter(default=1.5, significant=False,
                         description='Padding arount the solute (nm)')
    bt = luigi.Parameter(default="dodecahedron", significant=False,
                         description='Box type')
    salt_conc = luigi.FloatParameter(default=0.15, significant=False,
                         description='Salt (Na, Cl) concentration (mM)')
    n_repeats = luigi.IntParameter(default=3, significant=False,
                         description='Number of independent repeats')
    n_sampling_sims = luigi.IntParameter(default=1, significant=False,
                         description='Number of sampling simulations per independent repeat')
    
    def requires(self):
        return( Gather_Inputs_PL_folder(base_path=self.base_path,
                                        folder_path=self.folder_path,
                                        top_path=self.top_path,
                                        p=self.p, l=self.l) )
    
    def work(self):
        os.chdir(self.folder_path)
        
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
        
        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)
        
        #Return to basepath
        os.chdir(self.base_path)
        
    def output(self):
        files=[]
        #independent repeat (multiple for confidence estimate)
        for i in range(self.n_repeats):
            #sampling simulations in each repeat
            for m in range(self.n_sampling_sims):
                files.append("topol_ions%d_%d.top"%(i,m))
                files.append("ions%d_%d.pdb"%(i,m))
                
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


class Sim_EM(SGE_Sim):   
    #Signature by base_path, p, l, n, and i
    
    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    
    sim_path = luigi.Parameter(significant=False,
               description='Path where to run the simulation') #str
    
    def requires(self):
        return( Prep_PL_folder(base_path=self.base_path,
                               p=self.p, l=self.l) )
    
class Sim_NVT_posre(SGE_Sim):   
    #Signature by base_path, p, l, n, and i
    
    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    
    sim_path = luigi.Parameter(significant=False,
               description='Path where to run the simulation') #str
    
    def requires(self):
        return( Sim_EM(p=self.p, l=self.l,
                       i=self.i, m=self.m) )

# ==============================================================================
#                             Workflow Class
# ==============================================================================
class SGE_Workflow_alligned_inProtein(Workflow_alligned_inProtein):
                               
    def run_everything(self):
        """Runs the whole workflow.
        
        Parameters
        ----------
        None.
    
        Returns
        -------
        None.
        """
        
        pass
    
    