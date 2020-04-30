#!/usr/bin/env python

import luigi
import os
import shutil as sh
import subprocess
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.relFE.prep_folders_general import Gather_Inputs_folder_rel, Prep_folder_rel

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_PL_folder_rel(Gather_Inputs_folder):


    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.srctop="topol_abs_prot_norestr_amber.top"
        self.posre=True


class Prep_PL_folder_rel(Prep_folder): # will execute on the login node


    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_mdp="{}/protein/init.mdp".format(self.study_settings['mdp_path'])

    def requires(self):
        return([Gather_Inputs_PL_folder_rel(
                    folder_path=self.folder_path,
                    p=self.p, l=self.l,
                    study_settings=self.study_settings) ])
