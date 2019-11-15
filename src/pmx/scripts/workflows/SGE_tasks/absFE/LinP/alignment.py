#!/usr/bin/env python

import luigi
import numpy as np
import os
from luigi.contrib.sge import LocalSGEJobTask
from pmx import ndx
from pmx.model import Model
from pmx.scripts.workflows.fit_ligs_multiframes_python3 import fit,rotate_velocities_R
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.utils import read_from_mdp
from pmx.xtc import Trajectory


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_morphes(LocalSGEJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    sTI = luigi.Parameter(description='Coupling state for TI')

    folder_path = luigi.Parameter(significant=False,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    stage="morphes"

    #request 1 cores
    n_cpu = luigi.IntParameter(default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="", description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        description="Explicit job name given via qsub.")

    def __init__(self, *args, **kwargs):
        super(Task_PL_gen_morphes, self).__init__(*args, **kwargs)
        self._setupState()

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_test_{0}.mdp".format(
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
                  "> /dev/null 2>&1"%(tpr,trj,"frame.gro",self.study_settings['b']) )

        #restore base path
        os.chdir(self.base_path)


    def requires(self):
        return( Sim_PL_NPT(p=self.p, l=self.l, i=self.i, m=self.m, s=self.s,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )

    def output(self):
        #find nframes by reading the mdp file
        end_time, nstxout = read_from_mdp(self.mdp)
        nframes=int(end_time/nstxout) - int(self.study_settings['b']/nstxout)

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'frame%d.gro'%nf)) )
        return targets




class Task_PL_align(Task_PL_gen_morphes):

    def _setupState(self):
        if(not self.sTI == "C"):
            raise(ValueError("Aligning morphes for TI "
                     "state{}. ".format(self.sTI) +\
                     "Task_PL_align should only be done on TI stateC."))
            exit(1);
        else:
            self.s="B" # TI stateC depends on NPT sim in stateB


    def work(self):
        #make the C state
        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)

        m_A = Model(self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),bPDBTER=True)
        m_B = Model(self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),bPDBTER=True)
        m_C = Model(self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),bPDBTER=True)
        m_A.a2nm()
        m_B.a2nm()
        m_C.a2nm()

        srctraj=self.folder_path+"/state{2}/repeat{3}/npt{4}/traj.trr"

        trj_A = Trajectory(srctraj.format(self.p,self.l,"A",self.i,self.m))
        trj_B  = Trajectory(srctraj.format(self.p,self.l,"B",self.i,self.m))

        ndx_file = ndx.IndexFile(self.folder_path+"/index_prot_mol.ndx", verbose=False)
        p_ndx = np.asarray(ndx_file["Protein"].ids)-1
        l_ndx = np.asarray(ndx_file["MOL"].ids)-1

        #Frames are not acessible individually, just in sequence.
        #pmx.xtc.Trajectory is based on __iter__, so we need a custom
        #"for" loop to simultaneously go through both trajectories.
        #Based on https://www.programiz.com/python-programming/iterator
        iter_A = iter(trj_A)
        iter_B = iter(trj_B)
        fridx=0
        while True:
            try:
                frame_A = next(iter_A)
                frame_B = next(iter_B)
            except StopIteration:
                break

            if(not os.path.isfile("frame%d.gro"%fridx)):
                frame_A.update(m_A)
                frame_B.update(m_B)
                frame_B.update(m_C)

                # step1: fit prot from prot+lig onto apo protein
                (v1,v2,R) = fit( m_B, m_A, p_ndx, p_ndx )
                # rotate velocities
                # not needed. We aren't saving m_A

                # step2: ligand onto the ligand from prot+lig structure
                (v1,v2,R) = fit( m_A, m_C, l_ndx, l_ndx )
                # rotate velocities
                rotate_velocities_R( m_C, R )

                #replace coordinates and velocities of ligand in B with rotated ones from C
                for i in l_ndx:
                    m_B.atoms[i].x = m_C.atoms[i].x
                    m_B.atoms[i].v = m_C.atoms[i].v

                # output
                m_B.write("frame%d.gro"%fridx)

            fridx+=1

        #restore base path
        os.chdir(self.basepath)
