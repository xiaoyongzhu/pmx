#!/usr/bin/env python

import luigi
import os
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.restraints import Task_PL_gen_restraints
from pmx.scripts.workflows.find_anchors_and_write_ii_single_traj import find_restraints_align2crystal
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.alignment_align2crystal import Task_PL_align2crystal
from pmx.scripts.workflows.utils import check_file_ready


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_restraints_align2crystal(Task_PL_gen_restraints):

    #Parameters:
    restr_scheme = luigi.Parameter(significant=True,
                 description='Restraint scheme to use. '
                 'Aligned, Aligned_crystal, Fitted or Fixed')

    def __init__(self, *args, **kwargs):
        super(Task_PL_gen_restraints, self).__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']

        if(self.restr_scheme=="Aligned_crystal"):
            self.states=['D']
        else:
            raise(Exception(self.__class__.__name__+" only supports "
                            "the Aligned_crystal scheme."))

    def requires(self):
        reqs=[]
        #independent repeats for error analysis
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                reqs.append(Task_PL_align2crystal(
                                p=self.p, l=self.l, i=i, m=m, sTI='D',
                                study_settings=self.study_settings,
                                folder_path=self.folder_path,
                                parallel_env=self.parallel_env,
                                restr_scheme=self.restr_scheme))

        return(reqs)

    def output(self):
        targets=[
            luigi.LocalTarget(os.path.join(self.folder_path, 'index_prot_mol.ndx')),
            luigi.LocalTarget(os.path.join(self.folder_path, 'ii_aligned2crystal.itp')),
            luigi.LocalTarget(os.path.join(self.folder_path, 'out_dg_aligned2crystal.dat')),
            ]
        for s in self.states:
            targets.append([
                luigi.LocalTarget(os.path.join(self.folder_path, 'dump%s.gro'%s)),
                luigi.LocalTarget(os.path.join(self.folder_path, 'all_eq%s_fit.xtc'%s))
                ])
        for i in range(self.study_settings['n_repeats']):
            for m in range(self.study_settings['n_sampling_sims']):
                targets.append(luigi.LocalTarget(os.path.join(self.folder_path,
                                              "topolTI_aligned2crystal_ions%d_%d.top"%(i,m))))
        return targets

    def work(self):
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        srctpr=self.folder_path+"/state{2}/repeat{3}/npt{4}/tpr.tpr"
        srctraj=self.folder_path+"/state{2}/repeat{3}/aligned2crystal_morphes{4}/aligned.trr"

        #create prot+MOL index group
        os.system("echo \"1|13\nq\n\" | "
                  "gmx make_ndx -f ions0_0.pdb "
                  "-o index_prot_mol.ndx > /dev/null 2>&1")

        #make topology for prot+MOL
        os.system("sed 's/SOL/;SOL/g' topol.top > topol_prot_mol.top")

        #make topology for ApoP
        os.system("sed 's/SOL/;SOL/g' {base}/prot_{p}/apoP/topol.top "
                  "> topol_prot.top".format(base=self.base_path, p=self.p))

        s = 'D';

        #make tprs
        ref="box.pdb"
        ref_top="topol_prot_mol.top"
        mdp = self.study_settings['mdp_path'] + "/protein/init.mdp"

        os.system("gmx grompp -p {ref_top} -c {ref} -f {mdp} "
                  "-o tpr{s}.tpr -maxwarn 2 > grompp{s}.log 2>&1".format(
                      ref_top=ref_top, ref=ref, mdp=mdp, s=s) )
        check_file_ready("tpr{}.tpr".format(s))

        #collect trjs

        #remove previous log if it exists from a crashed attempt
        if(os.path.isfile("trjconv.log")):
            os.unlink("trjconv.log")

        #independent repeats for error analysis
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                tpr=srctpr.format(self.p,self.l,'A',i,m)
                trj=srctraj.format(self.p,self.l,s,i,m)
                sel="4 Protein_MOL"
                ndx="index_prot_mol.ndx"

                os.system("echo %s | "
                          "gmx trjconv -s %s -f %s "
                          "-o eq%s%d_%d.xtc "
                          "-sep -ur compact -pbc mol -center "
                          "-boxcenter zero -n %s "
                          "-b %d >> trjconv.log 2>&1"%(
                                  sel,tpr,trj, s,i,m, ndx,
                                  self.study_settings['b']) )

                check_file_ready("eq%s%d_%d.xtc"%(s,i,m))

        #concatenate trajectories
        os.system("gmx trjcat -f eq%s*.xtc -o all_eq%s.xtc -sort "
                  "-cat >> trjconv.log 2>&1"%(s,s) )
        check_file_ready("all_eq%s.xtc"%s)

        #fit to reference structure in tpr files
        os.system("echo 4 0 | gmx trjconv -s tpr%s.tpr -f all_eq%s.xtc "
                  "-o all_eq%s_fit.xtc -fit rot+trans "
                  ">> trjconv.log 2>&1"%(s,s,s) )
        check_file_ready("all_eq%s_fit.xtc"%s)

        #dump first frame
        os.system("echo 0 | gmx trjconv -f all_eq%s_fit.xtc "
                  "-s tpr%s.tpr -o dump%s.gro -dump 0 "
                  ">> trjconv.log 2>&1"%(s,s,s) )
        check_file_ready("dump%s.gro"%s)



        find_restraints_align2crystal(struct= 'dumpD.gro',
                        traj = "all_eqD_fit.xtc",
                        out="ii_aligned2crystal.itp",
                        an_cor_file="out_dg_aligned2crystal.dat",
                        log=False)


        #create a C state topology that holds ligand in place
        #with the restraint from the ii.itp files
        for i in range(self.study_settings['n_repeats']):
            for m in range(self.study_settings['n_sampling_sims']):
                top_ions="topol_ions%d_%d.top"%(i,m)
                topAC_ions="topolTI_aligned2crystal_ions%d_%d.top"%(i,m)
                topsrc=self.study_settings['top_path']+"/topol_abs_prot_restr_aligned2crystal_amber.top"
                if(not os.path.isfile(topsrc)):
                    raise(Exception("Missing reuired topology file: %s"%topsrc))
                os.system("cp {} {} > /dev/null 2>&1".format(top_ions,topAC_ions))
                os.system("tail -n 3 {} >> {}".format(topsrc,topAC_ions))
                check_file_ready(topAC_ions)

        #restore base path
        os.chdir(self.base_path)