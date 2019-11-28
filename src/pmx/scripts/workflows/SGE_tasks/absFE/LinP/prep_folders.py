#!/usr/bin/env python

import luigi
import os
import shutil as sh
import subprocess
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.absFE.prep_folders_general import Gather_Inputs_folder, Prep_folder
from pmx.scripts.workflows.SGE_tasks.absFE.ApoP.prep_folders import Prep_ApoP_folder

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_PL_folder(Gather_Inputs_folder):


    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.srctop="topol_abs_prot_norestr_amber.top"
        self.posre=True


class Prep_PL_folder(Prep_folder): # will execute on the login node


    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_mdp="{}/protein/init.mdp".format(self.study_settings['mdp_path'])

    def solvate(self):
        apoP_path=self.study_settings['base_path']+"/prot_{}/apoP".format(self.p)
        line = subprocess.check_output("tail -n 1 {}/topol_solvated.top".format(apoP_path),
                                       shell=True).decode('utf-8')
        if(not "SOL" in line):
            raise(Exception("Could not find the number of water molecules in ApoP."))
        need=int(line.split()[-1])

        d=self.study_settings['d']
        os.unlink("box.pdb")
        while(True):
            os.system("gmx editconf -f init.pdb -o box.pdb -bt %s -d %f "\
                      "> prep.log 2>&1"%(self.study_settings['bt'], d))

            os.system("gmx solvate -scale 1.0 -cp box.pdb -o water.pdb "\
                      "-cs spc216.gro -p topol_solvated.top -maxsol {} "
                      ">> prep.log 2>&1".format(need))

            line = subprocess.check_output("tail -n 1 topol_solvated.top",
                                           shell=True).decode('utf-8')
            have=int(line.split()[-1])

            if(have==need):
                break
            elif(have<need):
                #reset
                os.unlink("topol_solvated.top")
                os.unlink("box.pdb")
                os.unlink("water.pdb")
                sh.copy("topol.top","topol_solvated.top")
                d+=0.005 #add 0.05 A to the padding and try again
            else:
                raise(Exception("We have more water than expected. "
                                "This should never happen."))

    def requires(self):
        return([Gather_Inputs_PL_folder(
                    folder_path=self.folder_path,
                    p=self.p, l=self.l,
                    study_settings=self.study_settings),
                Prep_ApoP_folder(
                    folder_path=self.study_settings['base_path']+"/prot_{}/apoP".format(self.p),
                    p=self.p, study_settings=self.study_settings) ])
