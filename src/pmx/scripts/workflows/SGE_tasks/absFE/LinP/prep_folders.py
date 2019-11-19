#!/usr/bin/env python

import glob
import luigi
import os
import shutil as sh
from luigi.contrib.sge import LocalSGEJobTask
from pmx.scripts.workflows.utils import check_file_ready

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_PL_folder(LocalSGEJobTask): # will execute on the login node
    folder_path = luigi.Parameter(significant=False,
                      description='Path to the protein+ligand folder to set up')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

    study_settings = luigi.DictParameter(significant=False,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        description="Explicit job name given via qsub.")

    def work(self):

        #make folder
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        #topology
        sh.copy(self.study_settings['top_path']+"/topol_abs_prot_norestr_amber.top",
                self.folder_path+"/topol.top")
        sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/lig.itp",self.folder_path+"/lig.itp")
        sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.itp",self.folder_path+"/prot.itp")

        #initial coordinates where protein and ligand are bound
        sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_lig.pdb",
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
        os.chdir(self.study_settings['base_path'])

    def output(self):
        files=["topol.top", "lig.itp", "prot.itp", "init.pdb", "index.ndx",
               "prot_posre.itp", "prot_posre_soft.itp",
               "lig_posre.itp", "lig_posre_soft.itp"]
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


class Prep_PL_folder(LocalSGEJobTask): # will execute on the login node
    folder_path = luigi.Parameter(significant=False,
                         description='Path to the protein+ligand folder to set up')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

    study_settings = luigi.DictParameter(significant=False,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="", description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        description="Explicit job name given via qsub.")

    def requires(self):
        return( Gather_Inputs_PL_folder(folder_path=self.folder_path,
                                        p=self.p, l=self.l,
                                        study_settings=self.study_settings) )

    def work(self):
        os.chdir(self.folder_path)

        os.system("gmx editconf -f init.pdb -o box.pdb -bt %s -d %f "\
                  "> prep.log 2>&1"%(self.study_settings['bt'], self.study_settings['d']))
        check_file_ready("box.pdb")
        sh.copy("topol.top","topol_solvated.top")
        os.system("gmx solvate -scale 1.0 -cp box.pdb -o water.pdb "\
                  "-cs spc216.gro -p topol_solvated.top >> prep.log 2>&1")
        check_file_ready("water.pdb")
        os.system("gmx grompp -p topol_solvated.top -c water.pdb -o tpr.tpr "\
                  "-f %s/protein/init.mdp -v -maxwarn 2 "\
                  ">> prep.log 2>&1"%self.study_settings['mdp_path'])
        check_file_ready("tpr.tpr")

        #generate ions for each
        #independent repeat (multiple for confidence estimate)
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                top_ions="topol_ions%d_%d.top"%(i,m)
                pdb_ions="ions%d_%d.pdb"%(i,m)
                if(os.path.isfile(pdb_ions)): #skip if it already exists
                    continue
                sh.copy("topol_solvated.top", top_ions)
                os.system("echo 'SOL' | gmx genion -s tpr.tpr "
                          "-p %s -conc %f "
                          "-neutral -nname Cl -pname Na "
                          "-o %s >> genion.log 2>&1" %(
                              top_ions, self.study_settings['salt_conc'], pdb_ions) )
                check_file_ready(pdb_ions)

        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)

        #Return to basepath
        os.chdir(self.study_settings['base_path'])

    def output(self):
        files=[]
        #independent repeat (multiple for confidence estimate)
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                files.append("topol_ions%d_%d.top"%(i,m))
                files.append("ions%d_%d.pdb"%(i,m))

        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]
