#!/usr/bin/env python

import luigi
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.analysis import Task_PL_analysis_aligned
from pmx.scripts.workflows.SGE_tasks.relFE.LinW.TI import Task_PL_TI_simArray_rel


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_analysis_rel(Task_PL_analysis_aligned):

    #Parameters:
    folder_path = luigi.Parameter(significant=False,
        visibility=ParameterVisibility.HIDDEN,
        description='Path to the protein+ligand folder to set up')

    #request 1 core
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_l{l}_{i}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    
    restr_scheme = None

    def requires(self):
        tasks=[]

        for sTI in self.study_settings['TIstates']:
            for m in range(self.study_settings['n_sampling_sims']):
                tasks.append(Task_PL_TI_simArray_rel(
                      p=self.p, l=self.l, i=self.i, m=m, sTI=sTI,
                      study_settings=self.study_settings,
                      folder_path=self.folder_path,
                      parallel_env=self.parallel_env) )

        return(tasks)

