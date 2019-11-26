#!/usr/bin/env python

import glob
import luigi
import os
import pmx
import shutil as sh
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_folder(SGETunedJobTask):
    #run on the login node
    run_locally = luigi.BoolParameter(
        default = True,
        significant=False,
        parsing=luigi.BoolParameter.EXPLICIT_PARSING,
        description="run locally instead of on the cluster")

    #Parameters:
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

    #request 1 cores
    n_cpu = luigi.IntParameter(default=1, significant=False)
    parallel_env = luigi.Parameter(default='openmp_fast', significant=False)

    #variables to be overwriten in sub class' __init__()
    srctop=""
    posre=True

    def work(self):

        #make folder
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        #topology
        sh.copy(self.study_settings['top_path']+"/"+self.srctop, self.folder_path+"/topol.top")
        if(self.l):
            sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/lig.itp",self.folder_path+"/lig.itp")
        if(self.p):
            sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.itp",self.folder_path+"/prot.itp")

        #initial coordinates
        if(self.p and self.l): #P+L
            sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_"+self.l+".pdb",
                    self.folder_path+"/init.pdb")
        elif(self.p): #ApoP
            sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.pdb",
                    self.folder_path+"/init.pdb")
        elif(self.l): #L
            sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/ligand.pdb",
                    self.folder_path+"/init.pdb")

        #generate temporary index file
        os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                  "-o index.ndx > setup.log 2>&1")
        check_file_ready("index.ndx")

        if(self.posre):
            #generate restraints for equillibration
            #TODO: rewrite this to use the pmx Topology class
            if(self.p):
                os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                          "-fc 9000 9000 9000 -o prot_posre.itp "
                          "-n index.ndx >> setup.log 2>&1")
                check_file_ready("prot_posre.itp")
                os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                          "-fc 500 500 500 -o prot_posre_soft.itp "
                          "-n index.ndx >> setup.log 2>&1")
                check_file_ready("prot_posre_soft.itp")
            if(self.l):
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
        files=["topol.top", "init.pdb" , "index.ndx"]
        if(self.p):
            files.extend(["prot.itp"])
            if(self.posre):
                files.extend(["prot_posre.itp", "prot_posre_soft.itp"])
        if(self.l):
            files.extend(["lig.itp"])
            if(self.posre):
                files.extend(["lig_posre.itp", "lig_posre_soft.itp"])
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


class Prep_folder(SGETunedJobTask):
    #run on the login node
    run_locally = luigi.BoolParameter(
        default = True,
        significant=False,
        parsing=luigi.BoolParameter.EXPLICIT_PARSING,
        description="run locally instead of on the cluster")

    #Parameters:
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

    #request 1 cores
    n_cpu = luigi.IntParameter(default=1, significant=False)
    parallel_env = luigi.Parameter(default='openmp_fast', significant=False)

    def solvate(self):
        """Solvates the system.

        Subclasses can override this to set a custom number of water molecules

        Returns
        -------
        None.

        """
        os.system("gmx solvate -scale 1.0 -cp box.pdb -o water.pdb "\
                  "-cs spc216.gro -p topol_solvated.top >> prep.log 2>&1")


    def work(self):
        os.chdir(self.folder_path)

        os.system("gmx editconf -f init.pdb -o box.pdb -bt %s -d %f "\
                  "> prep.log 2>&1"%(self.study_settings['bt'], self.study_settings['d']))
        check_file_ready("box.pdb")
        sh.copy("topol.top","topol_solvated.top")

        #solvate using old VdW radii (found in mutff). This needs pmx's mutff
        #to be the first entry in GMXLIB. It contains the old version of
        #vdwradii.dat which produces better water density.
        if "GMXLIB" in os.environ.keys():
            orig_GMXLIB = os.environ["GMXLIB"]
        else:
            orig_GMXLIB = ""
        os.environ["GMXLIB"] = os.path.join(os.path.dirname(pmx.__file__),
                  "data/mutff/") + os.pathsep + orig_GMXLIB
        self.solvate()
        os.environ["GMXLIB"] = orig_GMXLIB

        check_file_ready("water.pdb")
        os.system("gmx grompp -p topol_solvated.top -c water.pdb -o tpr.tpr "\
                  "-f {} -v -maxwarn 2 "\
                  ">> prep.log 2>&1".format(self.init_mdp))
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
