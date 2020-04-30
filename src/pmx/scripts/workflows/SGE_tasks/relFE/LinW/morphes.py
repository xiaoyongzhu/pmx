#!/usr/bin/env python

import luigi
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.relFE.LinW.equil_sims import Sim_WL_NPT_rel
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.morphes import Task_PL_gen_morphes


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_WL_gen_morphes_rel(Task_PL_gen_morphes):

    #Parameters:
    p = None

    folder_path = luigi.Parameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Path to the water+ligand folder to set up')

    restr_scheme = luigi.Parameter(significant=False, default="",
                 description='Restraint scheme to use. '
                 'Aligned, Fitted or Fixed')

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #overwrite relevant variables
        self.mdp = self.study_settings['mdp_path'] +\
            "/water/eq_npt_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

    def requires(self):
        return( Sim_WL_NPT_rel(l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )


