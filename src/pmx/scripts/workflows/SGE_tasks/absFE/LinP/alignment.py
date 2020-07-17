#!/usr/bin/env python

import luigi
import numpy as np
import os
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx import ndx
from pmx.model import Model
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.fit_ligs_multiframes_python3 import fit,rotate_velocities_R, find_last_protein_atom
from pmx.scripts.workflows.utils import read_from_mdp
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.equil_sims import Sim_WL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.ApoP.equil_sims import Sim_ApoP_NPT
from pmx.xtc import Trajectory


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_align(SGETunedJobTask):
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

    extra_packages=[np]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._setupState()

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

        if(not (self.restr_scheme=="Aligned")):
            raise(Exception("{} should only be used if restraint "
                            "scheme is {}".format(
                                self.__class__.__name__ ,"Aligned")))

    def _setupState(self):
        if(self.sTI != "C"):
            raise(ValueError("Aligning morphes for TI state{}. "
                     "{} should only be done on TI stateC.".format(
                         self.sTI, self.__class__.__name__)))
            exit(1);
        # no need to set self.s as it is never used
        # else:
        #     self.s="B" # TI stateC depends on NPT sim in stateB

    def requires(self):
        #restraints require both state A & B for all repeats and sampling sims
        tasks=[]
        #for i in range(self.study_settings['n_repeats']):
            #for m in range(self.study_settings['n_sampling_sims']):

        #Currently work() is only susing self.i and self.m
        #TODO: all the m's should be concatenated and so align should depend
        #on all of them at once
        i=self.i
        m=self.m
        tasks.append(Sim_WL_NPT(l=self.l, i=i, m=m, s='B',
                  study_settings=self.study_settings,
                  folder_path=self.base_path+"/water/lig_{}".format(self.l),
                  parallel_env=self.parallel_env))
        tasks.append(Sim_PL_NPT(p=self.p, l=self.l, i=i, m=m, s='A',
                  study_settings=self.study_settings,
                  folder_path=self.folder_path,
                  parallel_env=self.parallel_env))
        tasks.append(Sim_ApoP_NPT(p=self.p, i=i, m=m,
                  study_settings=self.study_settings,
                  folder_path=self.base_path+"/prot_{}/apoP".format(self.p),
                  parallel_env=self.parallel_env))

        return(tasks)

    def output(self):
        #find nframes by reading the mdp file
        end_time, dtframe = read_from_mdp(self.mdp)
        nframes=int(end_time/dtframe) - int(self.study_settings['b']/dtframe) +1 #first frame counts

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'frame%d.gro'%nf)) )

        # targets.append(luigi.LocalTarget(
        #     os.path.join(self.sim_path, "aligned.trr")) )

        return targets

    def gen_ndx_w_chains(self, struct_path, ndx_path, nchains):
        #create default
        os.system("echo 'q' | gmx make_ndx -f {gro} -o {out} "
                  ">> align.log 2>&1".format(
                    gro=struct_path, out=ndx_path ) )
        #add chains
        n = ndx.IndexFile(ndx_path)
        m = Model(struct_path, bPDBTER=True)
        for i in range(nchains):
            ch = m.chains[i]
            grp = ndx.make_index_group(ch.atoms, "chain_" + ch.id)
            n.add_group(grp)
        n.write(ndx_path)


    def work(self):
    
        mylog=open("align.log","w")
    
        #find number of protein chains
        m_init = Model(self.folder_path+"/init.pdb", bPDBTER=True)
        n_prot_chains = len(m_init.chains)-1 #one is ligand
        if(n_prot_chains<1):
            raise(Exception("There should be at least one protein chain!"))

        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)
        if(os.path.isfile("align.log")): #clean old partial log if present
            os.unlink("align.log")

        trj_A_src=self.folder_path+"/state{2}/repeat{3}/npt{4}/".format(
            self.p, self.l, 'A', self.i, self.m)  #P+L
        trj_B_src=self.base_path+"/prot_{0}/apoP/repeat{3}/npt{4}/".format(
            self.p, self.l, None, self.i, self.m) #apoP
        trj_C_src=self.base_path+"/water/lig_{1}/state{2}/repeat{3}/npt{4}/".format(
            self.p, self.l, 'B', self.i, self.m)  #vacL

        #make the ndxs for the chains
        self.gen_ndx_w_chains(self.folder_path+"/ions%d_%d.pdb"%(self.i, self.m),
                              "PL_w_chains.ndx", n_prot_chains) #P+L
        check_file_ready("PL_w_chains.ndx")
        self.gen_ndx_w_chains(self.base_path+"/prot_{0}/apoP/ions{1}_{2}.pdb".format(
                                self.p, self.i, self.m),
                              "P_w_chains.ndx", n_prot_chains) #P
        check_file_ready("P_w_chains.ndx")
        #LW index
        os.system("echo 'q' | gmx make_ndx -f {gro} -o {out} "
                  ">> align.log 2>&1".format(
                  gro=self.base_path+"/water/lig_{1}/ions{3}_{4}.pdb".format(
                    self.p, self.l, 'B', self.i, self.m),
                  out="LW.ndx" ) )
        check_file_ready("LW.ndx")

        #Cut the begining off of trjs and center them
        names = ["A","B","C"]
        sources = [trj_A_src, trj_B_src, trj_C_src]
        ndxs = ["-n PL_w_chains.ndx", "-n P_w_chains.ndx", ""]
        for i in range(3):
            src=sources[i]

            #wrap mol centers
            os.system("echo System | gmx trjconv -s {tpr} -f {trj} -o {out} "
                  "-b {b} -ur compact -pbc mol"
                  ">> align.log 2>&1".format(
                      tpr=src+"tpr.tpr", trj=src+"traj.trr",
                      out="trj_{}_temp_cut_pbc.trr".format(names[i]),
                      b=self.study_settings['b']) )

            #center on chain_A and rewrap mol centers
            c = "chain_A"
            if(i==2): #LW
                c = "Other"
            os.system("echo {c} System | gmx trjconv -s {tpr} -f {trj} -o {out} "
                  "{ndx} -ur compact -center -pbc mol"
                  ">> align.log 2>&1".format(
                      tpr=src+"tpr.tpr", ndx=ndxs[i], c=c,
                      trj="trj_{}_temp_cut_pbc.trr".format(names[i]),
                      out="trj_{}.trr".format(names[i]) ) )

            check_file_ready("trj_{}.trr".format(names[i]))

            #clean temp
            os.unlink("trj_{}_temp_cut_pbc.trr".format(names[i]))

        #make the C state
        m_A = Model(self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),bPDBTER=True)
        m_B = Model(self.base_path+"/prot_{0}/apoP/ions{3}_{4}.pdb".format(
            self.p, self.l, None, self.i, self.m),bPDBTER=True) #apoP
        m_C = Model(self.base_path+"/water/lig_{1}/ions{3}_{4}.pdb".format(
            self.p, self.l, 'B', self.i, self.m),bPDBTER=True) #vacL
        m_A.a2nm()
        m_B.a2nm()
        m_C.a2nm()

        
        trj_A = Trajectory("trj_A.trr") #P+L
        trj_B = Trajectory("trj_B.trr") #apoP
        trj_C = Trajectory("trj_C.trr") #vacL



        ndx_file_A = ndx.IndexFile(self.folder_path+"/index_prot_mol.ndx", verbose=False)
        #ndx_file_C = ndx.IndexFile(self.base_path+"/water/lig_{}/index.ndx".format(self.l), verbose=False)
        ndx_file_C = ndx.IndexFile("LW.ndx", verbose=False)
        #p_ndx = np.asarray(ndx_file_A["Protein"].ids)-1
        p_ndx = np.asarray(ndx_file_A["C-alpha"].ids)-1 # as in Vytas' alignment script
        linA_ndx = np.asarray(ndx_file_A["MOL"].ids)-1
        l_ndx = np.asarray(ndx_file_C["MOL"].ids)-1
        
        
        #find chain and resID of the last residue of the protein
        mol_first_atom = m_A.atoms[linA_ndx[0]]
        chID = mol_first_atom.chain_id
        resID = mol_first_atom.resnr
        chain_local_res_index = -1;
        for i,r in enumerate(m_B.chdic[chID].residues):
            if(r.id==resID):
                chain_local_res_index=i;
                break;
        if(chain_local_res_index<0):
            raise("Could not find residue with resID %d in chain %s."%(last_prot_resID, chID))
        # chID,last_prot_resID = find_last_protein_atom( m_B )
        # #find the residue index to insert the ligand in the same chain as the end of the protein
        # chain_local_res_index = -1;
        # for i,r in enumerate(m_B.chdic[chID].residues):
            # if(r.id==last_prot_resID):
                # chain_local_res_index=i+1;
                # break;
        # if(chain_local_res_index<0):
            # raise("Could not find residue with resID %d in protein chain %s."%(last_prot_resID, chID))
        
        
        num_aligned_atoms = len(m_B.atoms) + l_ndx.shape[0]
        if(len(m_A.atoms) != num_aligned_atoms):
            print("\nWARNING: number of atoms in Apo + ligand ({}) does not match that of Holo ({})!\n".format(num_aligned_atoms,len(m_A.atoms)))
            mylog.write("\nWARNING: number of atoms in Apo + ligand ({}) does not match that of Holo ({})!\n".format(num_aligned_atoms,len(m_A.atoms)))
            mylog.flush()
        trj_out = Trajectory("aligned.trr", mode='Out',
                             atomNum = num_aligned_atoms) #aligned output

        #Frames are not acessible individually, just in sequence.
        #pmx.xtc.Trajectory is based on __iter__, so we need a custom
        #"for" loop to simultaneously go through both trajectories.
        #Based on https://www.programiz.com/python-programming/iterator
        iter_A = iter(trj_A)
        iter_B = iter(trj_B)
        iter_C = iter(trj_C)
        fridx=0
        mylog.write("\n\n\n\t\t###################\nStarting looping through frames:\n")
        mylog.flush()
        while True:
            try:
                frame_A = next(iter_A)
                frame_B = next(iter_B)
                frame_C = next(iter_C)
            except StopIteration:
                break

            #don't want frames while still equilibrating
            if(frame_A.time<self.study_settings['b']):
                continue

            if(not os.path.isfile("frame%d.gro"%fridx)):
                frame_A.update(m_A)
                #m_b needs to be reloaded to have correct # of atoms next iteration
                m_B = Model(self.base_path+"/prot_{0}/apoP/ions{3}_{4}.pdb".format(
                    self.p, self.l, None, self.i, self.m),bPDBTER=True) #apoP
                m_B.a2nm()
                frame_B.update(m_B, uv=True)
                frame_C.update(m_C, uv=True)
                mylog.write("frame %d: read models from pdbs\n"%fridx)
                #mylog.flush()

                # m_A.write("m_A1.gro")
                # step1: fit prot from prot+lig onto apo protein
                (v1,v2,R) = fit( m_B, m_A, p_ndx, p_ndx )
                # rotate velocities
                # not needed. We aren't saving m_A
                mylog.write("\t\tFit A on B\n")
                #mylog.flush()

                # step2: ligand onto the ligand from prot+lig structure
                (v1,v2,R) = fit( m_A, m_C, linA_ndx, l_ndx )
                mylog.write("\t\tFit C on A\n")
                #mylog.flush()
                # rotate velocities
                rotate_velocities_R( m_C, R )
                mylog.write("\t\tRotated C velocities\n")
                #mylog.flush()

                #insert vac ligand into B
                m_B.insert_residue(chain_local_res_index, m_C.residues[0], chID)
                mylog.write("\t\tOverwrote ligand in B with that from C\n")
                #mylog.flush()

                # #zero frame velocities so they don't get written to gro
                # for atom in m_B.atoms:
                #     for r in range(3):
                #         atom.v[r] = 0
                # output
                m_B.write("frame%d.gro"%fridx)
                mylog.write("\t\tWrote B as frame%d.gro\n"%fridx)
                mylog.flush()

                x = np.zeros(len(m_B.atoms)*3)
                v = np.zeros(len(m_B.atoms)*3)
                for i, atom in enumerate(m_B.atoms):
                    x[i*3:(i+1)*3]=atom.x
                    v[i*3:(i+1)*3]=atom.v

                # v=None
                mylog.write("\t\t\tSet atomx.x &.v; ready for writing trj frame to aligned.trr\n")
                trj_out.write_xtc_frame(step=frame_B.step, time=frame_B.time,
                                        lam=1.0, box=frame_B.box, x=x, v=v,
                                        units=m_B.unity, bTrr=True )

                mylog.write("\t\tWrote B to aligned trajectory\n")
                mylog.flush()


            fridx+=1

        trj_out.close()
        mylog.close()

        #restore base path
        os.chdir(self.base_path)
