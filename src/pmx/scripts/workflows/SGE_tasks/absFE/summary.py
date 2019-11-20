#!/usr/bin/env python

import copy
import luigi
import numpy as np
import matplotlib as plt
from luigi.contrib.sge import LocalSGEJobTask
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.analysis import Task_PL_analysis_aligned
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.analysis import Task_WL_analysis_aligned


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_summary_aligned(LocalSGEJobTask):

    #Parameters:
    hosts = luigi.ListParameter(description='list of protein names to evaluate')
    ligands = luigi.ListParameter(description='list of ligand names to evaluate')

    study_settings = luigi.DictParameter(significant=False,
        description='Dict of study stettings '
        'used to propagate settings to dependencies')

    #change default parallel environment
    parallel_env = luigi.Parameter(default='openmp_fast', significant=False)

    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        description="Explicit job name given via qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.WL_settings=copy.deepcopy(self.study_settings)
        self.WL_settings.update({'TIstates': self.WL_settings['states']})
        self.PL_settings=self.study_settings

    def work(self):
        #TODO: build summary table
        exit(1)
        pass

    def outputs(self):
        return([])

    def requires(self):
        tasks=[]

        #Ligand in Water
        p="water"
        for l in self.ligands:
            folder_path = self.basepath+'/'+p+'/lig_'+l
            for sTI in self.WL_settings['states']: #uses equil states for TI
                for i in range(self.WL_settings['n_repeats']):
                    tasks.append(Task_WL_analysis_aligned(
                        l = l, i = i,
                        study_settings = self.WL_settings,
                        folder_path = folder_path,
                        parallel_env=self.parallel_env))

        #Ligand in Protein
        for p in self.hosts:
            for l in self.ligands:
                folder_path = self.basepath+'/prot_'+p+'/lig_'+l
                for sTI in self.PL_settings['TIstates']:
                    for i in range(self.PL_settings['n_repeats']):
                        tasks.append(Task_PL_analysis_aligned(
                            p = p, l = l, i = i,
                            study_settings = self.PL_settings,
                            folder_path = folder_path,
                            parallel_env=self.parallel_env))

        return(tasks)

