#!/usr/bin/env python

import luigi
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.analysis import Task_PL_analysis_aligned
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.TI import Task_WL_TI_simArray


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_WL_analysis_aligned(Task_PL_analysis_aligned):

    #Parameters:
    p = None #disables base class' p

    folder_path = luigi.Parameter(significant=False,
        description='Path to the water+ligand folder to set up')

    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_l{l}_{i}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def requires(self):
        tasks=[]

        for sTI in self.study_settings['TIstates']:
            for m in range(self.study_settings['n_sampling_sims']):
                tasks.append(Task_WL_TI_simArray(
                      l=self.l, i=self.i, m=m, sTI=sTI,
                      study_settings=self.study_settings,
                      folder_path=self.folder_path,
                      parallel_env=self.parallel_env) )

        return(tasks)

