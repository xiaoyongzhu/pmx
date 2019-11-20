#!/usr/bin/env python

import luigi
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.equil_sims import Sim_WL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.alignment import Task_PL_gen_morphes


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_WL_gen_morphes(Task_PL_gen_morphes):

    #Parameters:
    p = luigi.Parameter(significant=False, default=None,
        description='Protein name') #disables base class' p

    folder_path = luigi.Parameter(significant=False,
                 description='Path to the water+ligand folder to set up')

    #request 1 cores
    n_cpu = luigi.IntParameter(default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #overwrite relevant variables
        self.mdp = self.study_settings['mdp_path'] +\
            "/water/eq_npt_test_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

    def requires(self):
        return( Sim_WL_NPT(l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )


