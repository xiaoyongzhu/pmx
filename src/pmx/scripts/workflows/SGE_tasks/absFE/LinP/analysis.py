#!/usr/bin/env python

import glob
import luigi
import os
import sys
import pmx.scripts.analyze_dhdl
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.TI import Task_PL_TI_simArray


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_analysis_aligned(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')

    folder_path = luigi.Parameter(significant=False,
        description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
        description='Dict of study stettings '
        'used to propagate settings to dependencies')

    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{i}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']
        self.ana_folder=self.folder_path+"/analysis/repeat%d"%self.i

    def work(self):
        #TODO: implement analysis
        exit(1)

        dHdlA=[]
        dHdlB=[]

        os.makedirs(self.ana_folder, exist_ok=True)
        os.chdir(self.ana_folder)

        for m in range(self.study_settings['n_sampling_sims']):
            self.folder_path+"/state%s/repeat%d/%s%d"%(
                        self.sTI, self.i, self.stage, self.m)
            dHdlA.append(glob.glob(self.folder_path +\
                        '/state{A}/repeat{i}/morphes{m}/dHdl*.xvg'.format(
                            A=self.study_settings['TIstates'][0],
                            i=self.i, m=m)))
            dHdlB.append(glob.glob(self.folder_path +\
                        '/state{B}/repeat{i}/morphes{m}/dHdl*.xvg'.format(
                            B=self.study_settings['TIstates'][1],
                            i=self.i, m=m)))

        #set script params and call analyze_dhdl
        orig_argv = sys.argv
        orig_stdout=sys.stdout
        orig_stderr=sys.stderr
        sys.argv = [['analyze_dhdl.py'],\
                    ['-fA'], dHdlA,
                    ['-fB'], dHdlB,
                    ['--nbins', "10"]]
        sys.argv = [item for sublist in sys.argv for item in sublist] #flatten argv

        with open("analysis.log", "w") as f:
            sys.stdout = f
            sys.stderr = f
            pmx.scripts.analyze_dhdl.entry_point()

        sys.argv = orig_argv #reset argv
        sys.stdout=orig_stdout #reset output
        sys.stderr=orig_stderr

        os.chdir(self.base_path)#reset cwd
        pass

    def outputs(self):
        files=['wplot.png', 'analysis.log']
        return([luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files])

    def requires(self):
        tasks=[]

        for sTI in self.study_settings['TIstates']:
            for m in range(self.study_settings['n_sampling_sims']):
                tasks.append(Task_PL_TI_simArray(
                      p=self.p, l=self.l, i=self.i, m=m, sTI=sTI,
                      study_settings=self.study_settings,
                      folder_path=self.folder_path,
                      parallel_env=self.parallel_env) )

        return(tasks)

