#!/usr/bin/env python

import glob
import luigi
import os
import shutil as sh
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.prep_folders import Gather_Inputs_PL_folder, Prep_PL_folder

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_WL_folder(Gather_Inputs_PL_folder):
    folder_path = luigi.Parameter(significant=False,
        description='Path to the water+ligand folder to set up')
    p = None #disables base class' p

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")



    def work(self):

        #make folder
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        #topology
        sh.copy(self.study_settings['top_path']+"/topol_abs_water_amber.top",
                self.folder_path+"/topol.top")
        sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/lig.itp",
                self.folder_path+"/lig.itp")

        #initial coordinates where protein and ligand are bound
        sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/ligand.pdb",
                self.folder_path+"/init.pdb")

        #generate temporary index file
        os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                  "-o index.ndx > setup.log 2>&1")
        check_file_ready("index.ndx")

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
        files=["topol.top", "lig.itp", "init.pdb", "index.ndx"]
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


class Prep_WL_folder(Prep_PL_folder): # will execute on the login node
    folder_path = luigi.Parameter(significant=False,
        description='Path to the water+ligand folder to set up')
    p = None #disables base class' p

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def requires(self):
        return( Gather_Inputs_WL_folder(folder_path=self.folder_path,
                                        l=self.l,
                                        study_settings=self.study_settings) )

    #work():    same as in PL
    #output():  same as in PL
