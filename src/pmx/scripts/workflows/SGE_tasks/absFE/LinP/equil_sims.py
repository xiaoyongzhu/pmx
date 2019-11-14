#!/usr/bin/env python

import luigi
import os
from pmx.scripts.workflows.SGE_tasks.Sim import SGE_Sim
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.prep_folders import Prep_PL_folder

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Sim_PL_EM(SGE_Sim):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    s = luigi.Parameter(description='Coupling state')

    folder_path = luigi.Parameter(significant=False,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(description='Dict of study stettings '
                      'used to propagate settings to dependencies')

    stage="em"
    #request 2 cores
    n_cpu = luigi.IntParameter(default=2, significant=False)

    def work(self):
        #set required file names

        states = self.study_settings['states']
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.s, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/em_posre_{0}.mdp".format(states[self.s])
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.base_path+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.base_path+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #call work in base class
        SGE_Sim.work(self)

    def requires(self):
        return( Prep_PL_folder(p=self.p, l=self.l,
                               study_settings=self.study_settings,
                               folder_path=self.folder_path) )
                                #no need to pass parallel_env as
                                #Prep_PL_folder runs on the login node
    def output(self):
        #output() is run before work()
        #so need to set sim_path here
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.s, self.i, self.stage, self.m)
        return luigi.LocalTarget(os.path.join(self.sim_path, 'confout.gro'))


class Sim_PL_NVT_posre(Sim_PL_EM):
    stage="nvt_posre"
    #request 4 cores
    n_cpu = luigi.IntParameter(default=4, significant=False)

    def work(self):
        #set required file names
        states = self.study_settings['states']
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.s, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_nvt_posre_{0}.mdp".format(states[self.s])
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.base_path+"/prot_{0}/lig_{1}/state{2}/repeat{3}/em{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.base_path+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #call work in base class
        SGE_Sim.work(self)

    def requires(self):
        return( Sim_PL_EM(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_PL_NVT_posre_soft(Sim_PL_EM):
    stage="nvt_posre_soft"
    #request 4 cores
    n_cpu = luigi.IntParameter(default=4, significant=False)

    def work(self):
        #set required file names
        states = self.study_settings['states']
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.s, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_nvt_posre_soft_{0}.mdp".format(states[self.s])
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.base_path+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.base_path+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #call work in base class
        SGE_Sim.work(self)

    def requires(self):
        return( Sim_PL_NVT_posre(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

class Sim_PL_NPT(Sim_PL_EM):
    stage="nvt_posre_soft"
    #request 4 cores
    n_cpu = luigi.IntParameter(default=4, significant=False)

    def work(self):
        #set required file names
        states = self.study_settings['states']
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.s, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_test_{0}.mdp".format(states[self.s])
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)
        self.struct = self.base_path+"/prot_{0}/lig_{1}/state{2}/repeat{3}/nvt_posre_soft{4}/confout.gro".format(
            self.p, self.l, self.s, self.i, self.m)
        self.posre = self.base_path+"/prot_{0}/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, self.s, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        #call work in base class
        SGE_Sim.work(self)

    def requires(self):
        return( Sim_PL_NVT_posre_soft(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )