#!/usr/bin/env python

import luigi
import numpy as np
import os
#from copy import deepcopy
from luigi.parameter import ParameterVisibility
from pmx import ndx
from pmx.model import Model
#from pmx.scripts.workflows.utils import read_from_mdp, bootstrap_Pearson, Pearson_err_func
from pmx.scripts.workflows.utils import read_from_mdp
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.restraints import Task_PL_gen_restraints
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.alignment import Task_PL_align
from pmx.xtc import Trajectory
from pmx.geometry import Rotation
#import pmx.scripts.workflows.postHoc_restraining_python3
#from pmx.scripts.workflows.postHoc_restraining_python3 import calc_dist, calc_angle, calc_dih, vector_prod, subtract_vecs, write_ii, write_dg
from pmx.scripts.workflows.postHoc_restraining_python3 import calc_dist, calc_angle, calc_dih, vector_prod, subtract_vecs
#from pmx.scripts.workflows.find_anchors_and_write_ii import wrap_ang, wrap_dih

#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#import scipy as sp


# ==============================================================================
#                               Helper Functions
# ==============================================================================
def readii(fii):
    lig=[-1,-1,-1]
    pro=[-1,-1,-1]
    ligfirst=True;
    nang=0
    ndih=0
    means=np.zeros(6)
    ks=np.zeros(6)
    with open(fii,'r') as f:
        block=""
        for cnt, line in enumerate(f):
            l=line.strip()
            if(not l): #empty line
                continue
            if('['in l): #read block name
                s=l.split()
                block=s[1]
                continue
            s=l.split()
            # if(block=="bonds"):
                # if(int(s[0])<int(s[1])):
                    # ligfirst=False

            #assume lig is first, we'll flip in the end if needed
            if(block=="bonds"):
                lig[0]=int(s[0])
                pro[0]=int(s[1])
                means[0]=float(s[3])
                ks[0]=float(s[-1])

            elif(block=="angles"):
                means[1+nang]=float(s[4])
                ks[1+nang]=float(s[-1])
                nang+=1

            elif(block=="dihedrals"):
                if(ndih==0):
                    lig=[int(i) for i in s[0:3]] #reverse order
                    lig=lig[::-1]
                elif(ndih==2):
                    pro=[int(i) for i in s[1:4]]

                means[3+ndih]=float(s[5])
                ks[3+ndih]=float(s[-1])
                ndih+=1

        # #flip lig & pro if not ligfirst
        #if(not ligfirst):
            #lig,pro=pro,lig

    return(np.array(lig), np.array(pro), means, ks) #1-indexed becasue bynum takes that

def update_anchors(model, idx_lig, idx_pro):
        L=[model.atoms[idx_lig[a]-1].x for a in range(3)]
        L.reverse()
        P=[model.atoms[idx_pro[a]-1].x for a in range(3)]
        return(L+P) #ligand first: [l2 l1 l0 p0 p2 p3]

def print_cur_rot(anchors, global_idcs, goal):
    order=[ #type, mean&sigma id, anchor indeces
            ["dih", 5, [5,4,3,2], "dih_C"], #dih_C [P2, P1, P0, L2]
            ["ang", 2,   [4,3,2], "ang_B"], #ang_B [P1, P0, P2]
            ["dist",0,     [3,2], "dist"],  #dist [P0, L2]
            ["dih", 4, [4,3,2,1], "dih_B"], #dih_B [P1, P0, L2, L1]
            ["ang", 1,   [3,2,1], "ang_A"], #ang_A [P0, L2, L1]
            ["dih", 3, [3,2,1,0], "dih_A"]  #dih_A [P0, L2, L1, L0]
          ]
    for op in order:
        cur_val=0
        indeces=op[2]
        if(op[0] == "dih"):
            cur_val = calc_dih(anchors,*indeces,False)
        elif(op[0] == "ang"):
            cur_val = calc_angle(anchors,*indeces,False)
        elif(op[0] == "dist"):
            cur_disp = subtract_vecs(anchors[indeces[1]],
                                     anchors[indeces[0]])
            cur_val = np.linalg.norm(cur_disp)

        print("\t{}:\t{:>7.4}  {:>7.4}\t{}".format(op[3], cur_val, goal[op[1]], [global_idcs[a] for a in indeces]))

def rotate_and_translate_Lig_to(g, lig, model, idx_lig, idx_pro, order=None):
    goal=g
    #goal[1]=np.pi-goal[1]
    #goal[2]=np.pi-goal[2]
    #goal[3]=-goal[3]
    #goal[4]=goal[4]
    #goal[5]=-goal[5]
    if(not order):
        order=[ #type, mean&sigma id, anchor indeces
                ["dih", 5, [5,4,3,2], "dih_C"], #dih_C [P2, P1, P0, L2]
                ["ang", 2,   [4,3,2], "ang_B"], #ang_B [P1, P0, P2]
                ["dist",0,     [3,2], "dist"],  #dist [P0, L2]
                ["dih", 4, [4,3,2,1], "dih_B"], #dih_B [P1, P0, L2, L1]
                ["ang", 1,   [3,2,1], "ang_A"], #ang_A [P0, L2, L1]
                ["dih", 3, [3,2,1,0], "dih_A"]  #dih_A [P0, L2, L1, L0]
              ]

    #print("lig:",idx_lig,"\tpro:",idx_pro)
    #global_idcs=idx_lig.tolist()[::-1]+idx_pro.tolist()
    #print("anchor indeces:",global_idcs)
    #print("goal:\n",goal)
    #print("original:")
    #anchors=update_anchors(model, idx_lig, idx_pro)
    #print_cur_rot(anchors, global_idcs, goal)

    opid=0
    for op in order:
        indeces=op[2]
        anchors=update_anchors(model, idx_lig, idx_pro)
        #print("op:", op)
        if(op[0] == "dih"):
            cur_val = calc_dih(anchors,*indeces,False)
            delta = goal[op[1]] - cur_val
            # delta = cur_val - goal[op[1]]
            # print("\tCurrent dih={}\t target dih={}\t diff={}".format(cur_val, goal[op[1]], delta))
            # print("\t\taxis points:",anchors[indeces[1]],anchors[indeces[2]])
            rot = Rotation(anchors[indeces[1]],anchors[indeces[2]])
            for num, a in enumerate(lig): #rot needs to be around 0,0,0, same as start of rot axis
                # print("\t\tidx=",num,"\tbefore rot:", a.x)
                a.x=rot.apply(a.x, delta)
                # print("\t\tidx=",num,"\tafter rot:", a.x)
                a.v=rot.apply(a.v, delta)
            # anchors=update_anchors(model, idx_lig, idx_pro)
            # print("\tdih after rotation={}\t target dih={}".format(calc_dih(anchors,*indeces,False), goal[op[1]]))
        elif(op[0] == "ang"):
            cur_val = calc_angle(anchors,*indeces,False)
            delta = goal[op[1]] - cur_val
            # print("\tCurrent ang={}\t target ang={}\t diff={}".format(cur_val, goal[op[1]], delta))
            normal = vector_prod(
                subtract_vecs(anchors[indeces[2]], anchors[indeces[1]]),
                subtract_vecs(anchors[indeces[1]], anchors[indeces[0]])
                )
            # print("\t\tnormal:",normal, "\tlen:", np.linalg.norm(normal))
            endp = np.array(anchors[indeces[1]])+np.array(normal)
            rot = Rotation(anchors[indeces[1]], endp)
            #for a in lig: #rot needs to be around 0,0,0, same as start of rot axis
            for num, a in enumerate(lig):
                #print("\t\tidx=",num,"\tbefore rot:", a.x)
                a.x=rot.apply(a.x, delta)
                #print("\t\tidx=",num,"\tafter rot:", a.x)
                a.v=rot.apply(a.v, delta)
            # anchors=update_anchors(model, idx_lig, idx_pro)
            # print("\tang after rotation={}\t target ang={}".format(calc_angle(anchors,*indeces,False), goal[op[1]]))
        elif(op[0] == "dist"):
            cur_disp = subtract_vecs(anchors[indeces[1]],
                                     anchors[indeces[0]])
            cur_val = np.linalg.norm(cur_disp)
            #print("\t\tcur_disp:",cur_disp, "\tlen:", cur_val)
            final_disp = cur_disp * goal[op[1]]/cur_val
            #print("\t\tfinal_disp:",final_disp, "\tlen:", np.linalg.norm(final_disp), "\tgoal:", goal[op[1]])
            delta = final_disp - cur_disp
            #print("\t\tdelta:", delta, "\tlen:", np.linalg.norm(delta))
            for a in lig:
                for ax in range(3):
                    a.x[ax] = a.x[ax] + delta[ax]
        else:
            raise(Exception("Unrecognized restraint type %s"%op[0]))

        #model.write("after_op_%d.gro"%opid)
        opid+=1

        #anchors=update_anchors(model, idx_lig, idx_pro)
        #print_cur_rot(anchors, global_idcs, goal)
        # exit(0);

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_decorelate_alg(SGETunedJobTask):
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

    T = luigi.FloatParameter(significant=False,
                 default=298.0, #K
                 description='Temperature in Kelvin. '
                 'Used to build distribution from restraint strengths.')

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

    #debug output
    debug = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default=False,
        description="show debug output in a log.")

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

        # if(not (self.restr_scheme=="Fitted")):
            # self.restr_scheme="Fitted"
            # raise(Exception("{} should only be used if restraint "
            #                 "scheme is {}".format(
            #                     self.__class__.__name__ ,"Fitted")))

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
        for i in range(self.study_settings['n_repeats']):
            for m in range(self.study_settings['n_sampling_sims']):
                tasks.append( Task_PL_gen_restraints(p=self.p, l=self.l,
                              i=self.i,
                              study_settings=self.study_settings,
                              folder_path=self.folder_path,
                              parallel_env=self.parallel_env,
                              restr_scheme=self.restr_scheme) )

                tasks.append(Task_PL_align(p=self.p, l=self.l, i=self.i, m=m, sTI='C',
                                  study_settings=self.study_settings,
                                  folder_path=self.folder_path,
                                  parallel_env=self.parallel_env,
                                  restr_scheme=self.restr_scheme) )

        return(tasks)

    def output(self):
        #find nframes by reading the mdp file
        end_time, dtframe = read_from_mdp(self.mdp)
        nframes=int(end_time/dtframe) - int(self.study_settings['b']/dtframe) +1 #first frame counts

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'start%d.gro'%nf)) )

        targets.append(luigi.LocalTarget(
            os.path.join(self.sim_path, "decorrelated.trr")) )

        # targets.append(luigi.LocalTarget(
        #         os.path.join(self.folder_path, "ii_{i}.itp".format(i=self.i))))
        # targets.append(luigi.LocalTarget(
        #         os.path.join(self.folder_path, "out_dg_{i}.dat".format(i=self.i))))

        return targets


    def work(self):
        #make the C state
        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)

        trj_A_src=self.sim_path+"/aligned.trr"  #P+L
        trj_A = Trajectory(trj_A_src) #P+L

        #m_A = Model(self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),bPDBTER=True)
        m_A = Model("frame0.gro",bPDBTER=False)
        m_A.a2nm()

        trj_out = Trajectory(self.sim_path+"/decorrelated.trr", mode='Out',
                             atomNum=len(m_A.atoms)) #aligned output

        ndx_file_A = ndx.IndexFile(self.folder_path+"/index_prot_mol.ndx", verbose=False)
        l_ndx = np.asarray(ndx_file_A["MOL"].ids)-1

        #select ligand
        lig=list(map(lambda i: m_A.atoms[i], l_ndx))

        #read restraints from ii and find decorelated sigmas
        # idx_lig, idx_pro, means, sigmas = self.find_decor_sigmas(
        #             self.folder_path+"/ii_cor_{i}.itp".format(i=self.i),
        #             self.folder_path+"/ions%d_%d.pdb"%(self.i,self.m),
        #             trj_A_src,
        #             plot=self.folder_path+"/correlations_{}.png".format(self.i) )
        # #means and sigmas in deg and nm



        idx_lig, idx_pro, means, ks = readii(self.folder_path+"/ii_{i}.itp".format(i=self.i))
        conv_ks=ks
        for k in range(1,6):
            conv_ks[k]=conv_ks[k]*(np.pi/180.0)**2
        kT = 8.31445985*0.001*self.T
        sigmas = np.sqrt(kT/conv_ks)
        #means and sigmas in deg and nm

        cov_mat=np.zeros((6,6))
        for j in range(6):
            if(j>0): #not the distance
                means[j]*=np.pi/180.0 #convert to rad
                sigmas[j]*=np.pi/180.0 #convert to rad
            cov_mat[j,j]=sigmas[j]*sigmas[j]

        if(self.debug):
            print("debug: starting on trajectory: {}".format(trj_A_src))

        #Frames are not acessible individually, just in sequence.
        #pmx.xtc.Trajectory is based on __iter__, so we need a custom
        #"for" loop to simultaneously go through both trajectories.
        #Based on https://www.programiz.com/python-programming/iterator
        iter_A = iter(trj_A)
        fridx=0
        while True:
            try:
                frame_A = next(iter_A)
            except StopIteration:
                break

            if(not os.path.isfile("start%d.gro"%fridx)):
                frame_A.update(m_A)
                if(self.debug):
                        print("debug: \tframe {}".format(fridx))
                #draw restraint coords from independent multivariate distribution
                sample = np.random.multivariate_normal(means, cov_mat) # nm & rad

                #wrap sample so that dihedrals are in (-180,180). Angles shouldn't need this.
                for k in range(3,len(sample)):
                    #if(k>=3): #dihedrals
                    #wrap sample so that dihedrals are in (-pi,pi). Angles shouldn't need this.
                    for k in range(3,len(sample)):
                        #if(k>=3): #dihedrals
                        if(sample[k]<-np.pi):
                           sample[k]+=2*np.pi
                        elif(sample[k]>np.pi):
                            sample[k]-=2*np.pi

                #rotate ligand to satisfy drawn restraint coords
                rotate_and_translate_Lig_to(sample, lig, m_A, idx_lig, idx_pro)

                # output
                x = np.zeros(len(m_A.atoms)*3)
                v = np.zeros(len(m_A.atoms)*3)
                for i, atom in enumerate(m_A.atoms):
                    x[i*3:(i+1)*3]=atom.x
                    v[i*3:(i+1)*3]=atom.v

                # v=None
                trj_out.write_xtc_frame(step=frame_A.step, time=frame_A.time,
                                        lam=1.0, box=frame_A.box, x=x, v=v,
                                        units=m_A.unity, bTrr=True )
                m_A.write("start%d.gro"%fridx)
            else:
                if(self.debug):
                    print("debug: \tframe {} exists".format(fridx))
            fridx+=1


        trj_out.close()
        if(self.debug):
            print("debug: done with trajectory".format(fridx))

        # decor_trjs=""
        # for m in range(self.study_settings['n_sampling_sims']):
        #     decor_trjs+=self.sim_path+"/frame*.gro"
        # ndxf="index_prot_mol_noH_{i}.ndx".format(i=self.i)
        # os.system("echo -e \"3\n22\n\" | python {script} "
        #           "-f {ap} "
        #           "-n {ndxf} -oii ii_{i}.itp -odg out_dg_{i}.dat > gen_restr_{i}.log 2>&1".format(
        #               ap=decor_trjs, ndxf=ndxf, i=self.i,
        #               script=pmx.scripts.workflows.postHoc_restraining_python3.__file__))

        # kT = 8.31445985*0.001*self.T
        # ks=kT/(np.array(sigmas)**2) #in kJ/mol * nm^-2 or rad^-2

        # #conver means back into deg
        # for j in range(1,6): #not the distance
        #     means[j]/=np.pi/180.0 #convert to deg

        # #write new ii
        # ii = {}
        # lig1 = idx_lig[0]
        # lig2 = idx_lig[1]
        # lig3 = idx_lig[2]
        # prot1 = idx_pro[0]
        # prot2 = idx_pro[1]
        # prot3 = idx_pro[2]

        # bond1 = [ prot1, lig1, 6, [means[0], 0.0, means[0], ks[0]] ]
        # ii['bonds'] = [bond1]
        # angle1 = [ prot2, prot1, lig1, 1, [means[1], 0.0, means[1],  ks[1]] ]
        # angle2 = [ prot1, lig1, lig2, 1, [means[2], 0.0, means[2],  ks[2]] ]
        # ii['angles'] = [ angle1, angle2 ]
        # dihedral1 = [ prot3, prot2, prot1, lig1, 2, [means[3], 0.0, means[3],  ks[3]] ]
        # dihedral2 = [ prot2, prot1, lig1, lig2, 2, [means[4], 0.0, means[4],  ks[4]] ]
        # dihedral3 = [ prot1, lig1, lig2, lig3, 2, [means[5], 0.0, means[5],  ks[5]] ]
        # ii['dihedrals'] = [ dihedral1, dihedral2, dihedral3 ]
        # write_ii( ii, self.folder_path+"/ii_{i}.itp".format(i=self.i) )

        # # calculate dG contribution
        # V0 = 1.66            # standard volume in nm^3
        # dgPrefactor = ( 8.0*np.power(np.pi,2.0)*V0/(np.power(means[0],2.0)*np.sin(means[1]*np.pi/180.0)*np.sin(means[2]*np.pi/180.0)) )
        # dgForceConstants = np.sqrt(ks[0]*ks[1]*ks[2]*ks[3]*ks[4]*ks[5])/np.power(2.0*np.pi*kT,3.0)
        # dg = -kT * np.log(dgPrefactor*dgForceConstants)
        # dg = np.round(dg,2)

        # # output dG contribution
        # write_dg( dg, self.folder_path+"/out_dg_{i}.dat".format(i=self.i) )


        #restore base path
        os.chdir(self.base_path)



    # def find_decor_sigmas(self, ii, struct, traj, plot=None):
    #     idx_lig, idx_pro, means, ks = readii(ii)

    #     m = Model(struct)
    #     m.a2nm() #set units to nm
    #     trj_m = Trajectory(traj) #P+L
    #     iter_m = iter(trj_m)

    #     fr=0
    #     dist=[]
    #     angA=[]
    #     angB=[]
    #     dihA=[]
    #     dihB=[]
    #     dihC=[]
    #     while True:
    #         try:
    #             frame_m = next(iter_m)
    #         except StopIteration:
    #             break

    #         frame_m.update(m)
    #         anchors=update_anchors(m, idx_lig, idx_pro)

    #         dist.append(calc_dist(  anchors, 3, 2 )) #nm
    #         angA.append(calc_angle( anchors, 4, 3, 2)) #deg
    #         angB.append(calc_angle( anchors, 3, 2, 1)) #deg

    #         dihA.append(calc_dih(   anchors, 5, 4, 3, 2)) #deg
    #         dihB.append(calc_dih(   anchors, 4, 3, 2, 1)) #deg
    #         dihC.append(calc_dih(   anchors, 3, 2, 1, 0)) #deg

    #         fr+=1
    #     # print("went through {} frames for repeat {}".format(fr, self.i))

    #     #make a coarse histogram to find the peak
    #     bscale=10 # deg/bin
    #     angA_h = np.histogram(angA, bins=int(180/bscale), range=(0,180), density=True)
    #     angB_h = np.histogram(angB, bins=int(180/bscale), range=(0,180), density=True)
    #     dihA_h = np.histogram(dihA, bins=int(360/bscale), range=(-180,180), density=True)
    #     dihB_h = np.histogram(dihB, bins=int(360/bscale), range=(-180,180), density=True)
    #     dihC_h = np.histogram(dihC, bins=int(360/bscale), range=(-180,180), density=True)


    #     dists=[None, angA_h, angB_h, dihA_h, dihB_h, dihC_h]
    #     data=[dist, angA, angB, dihA, dihB, dihC]
    #     for i in range(len(data)):
    #         data[i]=np.array(data[i], dtype=np.float64)

    #         #wrap values to center them on the peakdists=[dist_h, angA_h, angB_h, dihA_h, dihB_h, dihC_h]
    #         if(i>0):
    #             d=dists[i]
    #             dists[i] = [d[0], d[1][:-1]+(d[1][1]-d[1][0])*0.5]
    #             m=dists[i][1][np.argmax(dists[i][0])] #coord of max height of distribution
    #             if(i<3): #wrap angles
    #                 data[i]=wrap_ang(data[i],m)
    #             else:    #wrap dihedrals
    #                 data[i]=wrap_dih(data[i],m)

    #     data[0]*=100.0 #nm * 10^-2 to avoid floating point issues in eigh()

    #     new_sigmas=np.std(data, axis=-1)
    #     new_means=np.mean(data, axis=-1)

    #     new_vars=new_sigmas*new_sigmas
    #     for i in range(6):
    #         for j in range(i+1,6):
    #             values=np.vstack([data[i],data[j]])
    #             cov=np.cov(values)
    #             w, v = np.linalg.eigh(cov)
    #             xvar = 1./np.sqrt(v[1,0]**2/w[1]**2 + v[1,1]**2/w[0]**2)
    #             yvar = 1./np.sqrt(v[0,0]**2/w[1]**2 + v[0,1]**2/w[0]**2)
    #             xsig = np.sqrt(xvar)
    #             ysig = np.sqrt(yvar)
    #             if(xsig<new_sigmas[i]):
    #                 new_sigmas[i]=xsig
    #                 new_vars[i]=xvar
    #             if(ysig<new_sigmas[j]):
    #                 new_sigmas[j]=ysig
    #                 new_vars[j]=yvar

    #     dist_names=["distance", "angA", "angB", "dihA", "dihB", "dihC"]
    #     units=["10^-2 nm", "deg", "deg", "deg", "deg", "deg"]
    #     # units_FC=["kJ/(mol*nm^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)"]

    #     if(plot):
    #         plt.figure(figsize=(18, 18))
    #         sp_grid = gridspec.GridSpec(5, 5)

    #     for i in range(6):
    #         for j in range(i+1,6):
    #             xspread=np.max(data[i])-np.min(data[i])
    #             xlim=(np.min(data[i])-0.2*xspread, np.max(data[i])+0.2*xspread)
    #             yspread=np.max(data[j])-np.min(data[j])
    #             ylim=(np.min(data[j])-0.2*yspread, np.max(data[j])+0.2*yspread)

    #             if(plot):
    #                 #Pearson correlation coef
    #                 pearson = Pearson_err_func(data[i],data[j])
    #                 pearson_err = bootstrap_Pearson(data[i],data[j])

    #                 ax=plt.subplot(sp_grid[j-1,i])

    #                 ax.set_xlim(xlim)
    #                 ax.set_ylim(ylim)

    #                 plt.title("%s vs %s"%(dist_names[i],dist_names[j]))
    #                 plt.xlabel("%s (%s)"%(dist_names[i],units[i]))
    #                 plt.ylabel("%s (%s)"%(dist_names[j],units[j]))

    #                 #points
    #                 plt.plot(data[i], data[j], '.k')

    #                 #new pdf
    #                 new_cov=np.zeros((2,2))
    #                 new_cov[np.diag_indices_from(new_cov)] = [new_vars[i],new_vars[j]]
    #                 new_distr = sp.stats.multivariate_normal([new_means[i],new_means[j]], new_cov)
    #                 X, Y = np.mgrid[xlim[0]:xlim[1]:100j,
    #                                 ylim[0]:ylim[1]:100j]
    #                 pos = np.empty(X.shape + (2,))
    #                 pos[:, :, 0] = X; pos[:, :, 1] = Y
    #                 new_Z = new_distr.pdf(pos)
    #                 ax.imshow(np.rot90(new_Z), cmap=plt.cm.gist_earth_r,
    #                          extent=[xlim[0],xlim[1],ylim[0],ylim[1]],
    #                          aspect='auto')

    #                 mean=[np.mean(data[i]), np.mean(data[j])]
    #                 ax.plot(mean[0], mean[1], 'r+', markersize=10)

    #                 # ax.arrow(*mean, *new_cov[1], color='r', width=0.08)
    #                 # ax.arrow(*mean, *new_cov[0], color='b', width=0.08)

    #                 ax.annotate("Pearson CC=%.3f+-%.3f"%(pearson, pearson_err),
    #                         xy=(0.98, 0.02),
    #                         xycoords='axes fraction', textcoords='offset points',
    #                         size='large', ha='right', va='baseline')

    #     if plot:
    #         plt.tight_layout()
    #         plt.savefig(plot, dpi=300)

    #     #convert distance back into nm and return distribution widths
    #     new_sigmas[0]/=100.
    #     new_means[0]/=100.

    #     return(idx_lig, idx_pro, new_means, new_sigmas)

