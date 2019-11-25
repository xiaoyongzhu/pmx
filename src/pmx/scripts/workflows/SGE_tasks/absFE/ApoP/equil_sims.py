#!/usr/bin/env python

import luigi
from pmx.scripts.workflows.SGE_tasks.absFE.ApoP.prep_folders import Prep_ApoP_folder
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_EM, Sim_PL_NVT_posre, Sim_PL_NVT_posre_soft, Sim_PL_NPT

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Sim_ApoP_EM(Sim_PL_EM):

    #Parameters:
    l = None #disables base class' l
    s = None

    folder_path = luigi.Parameter(significant=False,
                 description='Path to the protein+ligand folder to set up')

    #request 2 cores
    n_cpu = luigi.IntParameter(default=2, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super(Sim_PL_EM, self).__init__(*args, **kwargs)
        #set required file names
        self.base_path = self.study_settings['base_path']
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/apo_protein/em_posre.mdp"
        self.sim_path = self.folder_path+"/repeat{i}/{stage}{m}".format(
            i=self.i, stage=self.stage, m=self.m)

    def requires(self):
        return( Prep_ApoP_folder(p=self.p,
                               study_settings=self.study_settings,
                               folder_path=self.folder_path) )
                                #no need to pass parallel_env as
                                #Prep_WL_folder runs on the login node


class Sim_ApoP_NVT_posre(Sim_PL_NVT_posre):
    #Parameters:
    l = None #disables base class' l
    s = None

    #request 4 cores
    n_cpu = luigi.IntParameter(default=4, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super(Sim_PL_EM, self).__init__(*args, **kwargs)
        #set required file names
        self.base_path = self.study_settings['base_path']
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/apo_protein/eq_nvt_posre.mdp"
        self.struct = self.folder_path+"/repeat{i}/em{m}/confout.gro".format(
            i=self.i, m=self.m)
        self.sim_path = self.folder_path+"/repeat{i}/{stage}{m}".format(
            i=self.i, stage=self.stage, m=self.m)

    def requires(self):
        return( Sim_ApoP_EM(p=self.p, i=self.i, m=self.m,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_ApoP_NVT_posre_soft(Sim_PL_NVT_posre_soft):
    #Parameters:
    l = None #disables base class' l
    s = None

    #request 4 cores
    n_cpu = luigi.IntParameter(default=4, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super(Sim_PL_EM, self).__init__(*args, **kwargs)
        #set required file names
        self.base_path = self.study_settings['base_path']
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/apo_protein/eq_nvt_posre_soft.mdp"
        self.struct = self.folder_path+"/repeat{i}/nvt_posre{m}/confout.gro".format(
            i=self.i, m=self.m)
        self.sim_path = self.folder_path+"/repeat{i}/{stage}{m}".format(
            i=self.i, stage=self.stage, m=self.m)

    def requires(self):
        return( Sim_ApoP_NVT_posre(p=self.p, i=self.i, m=self.m,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )


class Sim_ApoP_NPT(Sim_PL_NPT):
    #Parameters:
    l = None #disables base class' l
    s = None

    #request 4 cores
    n_cpu = luigi.IntParameter(default=4, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super(Sim_PL_EM, self).__init__(*args, **kwargs)
        #set required file names
        self.base_path = self.study_settings['base_path']
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.folder_path+"/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #override relevant file names
        self.mdp = self.study_settings['mdp_path'] +\
            "/apo_protein/eq_npt_test.mdp"
        self.struct = self.folder_path+"/repeat{i}/nvt_posre_soft{m}/confout.gro".format(
            i=self.i, m=self.m)
        self.sim_path = self.folder_path+"/repeat{i}/{stage}{m}".format(
            i=self.i, stage=self.stage, m=self.m)



    def requires(self):
        return( Sim_ApoP_NVT_posre_soft(p=self.p, i=self.i, m=self.m,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )