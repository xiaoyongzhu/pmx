#!/usr/bin/env python

from pmx.scripts.workflows.SGE_tasks.absFE.LinP.alignment import Task_PL_gen_morphes
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_morphes_align2crystal(Task_PL_gen_morphes):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #overwrite variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/aligned2crystal_%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

    def requires(self):
        return( Sim_PL_NPT(p=self.p, l=self.l, i=self.i, m=self.m, s='A',
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )
