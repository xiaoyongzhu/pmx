#!/usr/bin/env python

import luigi
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.absFE.prep_folders_general import Gather_Inputs_folder, Prep_folder

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_WL_folder(Gather_Inputs_folder):
    p = None #disables base class' p

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.srctop="topol_abs_water_amber.top"
        self.posre=False


class Prep_WL_folder(Prep_folder): # will execute on the login node
    p = None #disables base class' p

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.init_mdp="{}/water/init.mdp".format(self.study_settings['mdp_path'])

    def requires(self):
        return( Gather_Inputs_WL_folder(folder_path=self.folder_path,
                                        l=self.l,
                                        study_settings=self.study_settings) )

