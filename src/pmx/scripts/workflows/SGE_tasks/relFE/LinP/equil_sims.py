#!/usr/bin/env python

import luigi
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.relFE.LinP.prep_folders import Prep_PL_folder_rel
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_EM, Sim_PL_NVT_posre, Sim_PL_NPT

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Sim_PL_EM_rel(Sim_PL_EM):

    def requires(self):
        return( Prep_PL_folder_rel(p=self.p, l=self.l,
                               study_settings=self.study_settings,
                               folder_path=self.folder_path) )
                                #no need to pass parallel_env as
                                #Prep_PL_folder runs on the login node


class Sim_PL_NVT_posre_rel(Sim_PL_NVT_posre):

    def requires(self):
        return( Sim_PL_EM_rel(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_PL_NVT_posre_soft_rel(Sim_PL_NVT_posre_soft):

    def requires(self):
        return( Sim_PL_NVT_posre_rel(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_PL_NPT_rel(Sim_PL_NPT):
    
    def requires(self):
        return( Sim_PL_NVT_posre_soft_rel(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )
