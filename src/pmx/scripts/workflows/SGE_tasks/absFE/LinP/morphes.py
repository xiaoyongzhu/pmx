#!/usr/bin/env python

import glob
import luigi
import os
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.restraints import Task_PL_gen_restraints
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.utils import read_from_mdp


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_morphes(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    sTI = luigi.Parameter(description='Coupling state for TI')

    folder_path = luigi.Parameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    restr_scheme = luigi.Parameter(significant=True,
                 description='Restraint scheme to use. '
                 'Aligned, Fitted or Fixed')

    stage="morphes"

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super(Task_PL_gen_morphes, self).__init__(*args, **kwargs)
        self._setupState()

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

    def _setupState(self):
        self.s=self.sTI

    def work(self):
        #generate morphs for A state
        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)

        tpr=self.folder_path+"/state{2}/repeat{3}/npt{4}/tpr.tpr".format(
            self.p, self.l, self.s, self.i, self.m)
        trj=self.folder_path+"/state{2}/repeat{3}/npt{4}/traj.trr".format(
            self.p, self.l, self.s, self.i, self.m)

        #this is slow
        os.system("echo 0 | gmx trjconv -s %s "
                  "-f %s -o %s "
                  "-b %f -sep -ur compact -pbc mol "
                  "> trjconv.log 2>&1"%(tpr,trj,"frame.gro",self.study_settings['b']) )

        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)

        #restore base path
        os.chdir(self.base_path)


    def requires(self):
        #restraints require both state A & B for all repeats and sampling sims
        # return( Task_PL_gen_restraints(p=self.p, l=self.l,
                          # i=self.i,
                          # study_settings=self.study_settings,
                          # folder_path=self.folder_path,
                          # parallel_env=self.parallel_env,
                          # restr_scheme=self.restr_scheme) )
                          
        return( Sim_PL_NPT(p=self.p, l=self.l, i=self.i, m=self.m, s=sTI,
                                  study_settings=self.study_settings,
                                  folder_path=self.folder_path,
                                  parallel_env=self.parallel_env) )

    def output(self):
        #find nframes by reading the mdp file
        end_time, dtframe = read_from_mdp(self.mdp)
        nframes=int(end_time/dtframe) - int(self.study_settings['b']/dtframe) +1 #first frame counts

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'frame%d.gro'%nf)) )
        return targets

