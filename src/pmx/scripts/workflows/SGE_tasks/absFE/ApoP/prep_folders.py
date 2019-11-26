#!/usr/bin/env python

import luigi
from pmx.scripts.workflows.SGE_tasks.absFE.prep_folders_general import Gather_Inputs_folder, Prep_folder

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_ApoP_folder(Gather_Inputs_folder):
    l = None #disables base class' l

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.srctop="topol_abs_apo_protein_amber.top"
        self.posre=True


class Prep_ApoP_folder(Prep_folder): # will execute on the login node
    l = None #disables base class' l

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_mdp="{}/apo_protein/init.mdp".format(self.study_settings['mdp_path'])

    def requires(self):
        return( Gather_Inputs_ApoP_folder(folder_path=self.folder_path,
                                        p=self.p,
                                        study_settings=self.study_settings) )

