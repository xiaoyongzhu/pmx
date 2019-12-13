#!/usr/bin/env python

import luigi
import MDAnalysis as md
import os
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.find_anchors_and_write_ii import find_restraints
from pmx.scripts.workflows.find_anchors_and_write_ii_single_traj import find_restraints_align2crystal
from pmx.scripts.workflows.find_avg import find_avg_struct
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.ApoP.equil_sims import Sim_ApoP_NPT
from pmx.scripts.workflows.utils import check_file_ready


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_restraints(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

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

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    extra_packages=[md]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']

        if(self.restr_scheme=="Aligned" or self.restr_scheme=="Fitted"):
            self.states=["A", "ApoProt"]
        elif(self.restr_scheme=="Fixed"):
            raise(Exception("Fixed restraints not yet implemened."))
            self.states = self.study_settings['states']
        else:
            raise(Exception("Unrecognized restraint scheme '%s'"%self.restr_scheme))

    def work(self):
        #generate morphs for A state
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        srctpr=self.folder_path+"/state{2}/repeat{3}/npt{4}/tpr.tpr"
        srctraj=self.folder_path+"/state{2}/repeat{3}/npt{4}/traj.trr"

        # #create prot+MOL index group
        # os.system("echo \"1|13\nq\n\" | "
        #           "gmx make_ndx -f ions0_0.pdb "
        #           "-o index_prot_mol.ndx > /dev/null 2>&1")

        #make topology for prot+MOL
        os.system("sed 's/SOL/;SOL/g' topol.top > topol_prot_mol.top")

        #make topology for ApoP
        os.system("sed 's/SOL/;SOL/g' {base}/prot_{p}/apoP/topol.top "
                  "> topol_prot.top".format(base=self.base_path, p=self.p))

        for s in self.states:
            #make tprs
            if(s == "A"):   #align A to initial structure
                ref="box.pdb"
                ref_top="topol_prot_mol.top"
                mdp = self.study_settings['mdp_path'] + "/protein/init.mdp"
            else:           #align B to average of A
                ref="averageA_prot_only.gro"
                ref_top="topol_prot.top"
                mdp = self.study_settings['mdp_path'] + "/apo_protein/init.mdp"


            os.system("gmx grompp -p {ref_top} -c {ref} -f {mdp} "
                      "-o tpr{s}.tpr -maxwarn 2 > grompp{s}.log 2>&1".format(
                          ref_top=ref_top, ref=ref, mdp=mdp, s=s) )
            check_file_ready("tpr{}.tpr".format(s))

            #collect trjs
            # print("\tCollecting trajectories for state%s"%s)

            #remove previous log if it exists from a crashed attempt
            if(os.path.isfile("trjconv.log")):
                os.unlink("trjconv.log")

            #independent repeats for error analysis
            for i in range(self.study_settings['n_repeats']):
                #sampling simulations in each repeat
                for m in range(self.study_settings['n_sampling_sims']):
                    if(s=="ApoProt"):
                        apoP_path=self.base_path+"/prot_{0}/apoP/repeat{3}/npt{4}/"
                        tpr=apoP_path.format(self.p,self.l,s,i,m)+"tpr.tpr"
                        trj=apoP_path.format(self.p,self.l,s,i,m)+"traj.trr"
                        sel="4 Protein"
                        ndx=self.base_path+"/prot_{0}/apoP/index.ndx".format(self.p)
                    else:
                        tpr=srctpr.format(self.p,self.l,s,i,m)
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

            #find avg structure of A
            if(s=="A"):
                # print("\tFinding average structure")
                find_avg_struct("dumpA.gro", "all_eqA_fit.xtc",
                                "averageA.gro")
                check_file_ready("averageA.gro")
                # print("\tExtracting prot. only")

                os.system("echo Protein | gmx trjconv -s tprA.tpr -f averageA.gro "
                          "-o averageA_prot_only.gro >> trjconv.log 2>&1" )
                check_file_ready("averageA_prot_only.gro")


        #generate the restraints
        # print("\tGenerating the restraints")
        if(self.restr_scheme=="Aligned" or self.restr_scheme=="Fitted"):
            find_restraints(structB= "dumpApoProt.gro",
                            trajB = "all_eqApoProt_fit.xtc",
                            log=False)
        elif(self.restr_scheme=="Fitted"):
            find_restraints(log=False)

        check_file_ready(os.path.join('ii.itp'))

        #create a C state topology that holds ligand in place
        #with the restraint from the ii.itp files
        for i in range(self.study_settings['n_repeats']):
            for m in range(self.study_settings['n_sampling_sims']):
                top_ions="topol_ions%d_%d.top"%(i,m)
                topAC_ions="topolTI_ions%d_%d.top"%(i,m)
                topsrc=self.study_settings['top_path']+"/topol_abs_prot_restr_amber.top"
                os.system("cp {} {} > /dev/null 2>&1".format(top_ions,topAC_ions))
                os.system("tail -n 3 {} >> {}".format(topsrc,topAC_ions))
                check_file_ready(topAC_ions)

        #restore base path
        os.chdir(self.base_path)


    def requires(self):
        reqs=[]
        #independent repeats for error analysis
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                #states of equilibrium sims
                if(self.restr_scheme=="Aligned" or self.restr_scheme=="Fitted"):
                    reqs.append(Sim_PL_NPT(p=self.p, l=self.l, i=i, m=m, s='A',
                                      study_settings=self.study_settings,
                                      folder_path=self.folder_path,
                                      parallel_env=self.parallel_env)
                            )
                    reqs.append(Sim_ApoP_NPT(p=self.p, i=i, m=m,
                          study_settings=self.study_settings,
                          folder_path=self.base_path+"/prot_{}/apoP".format(self.p),
                          parallel_env=self.parallel_env))
                else:
                    for s in self.states:
                        reqs.append(
                            Sim_PL_NPT(p=self.p, l=self.l, i=i, m=m, s=s,
                                      study_settings=self.study_settings,
                                      folder_path=self.folder_path,
                                      parallel_env=self.parallel_env)
                            )


        return(reqs)

    def output(self):
        targets=[
            luigi.LocalTarget(os.path.join(self.folder_path, 'index_prot_mol.ndx')),
            luigi.LocalTarget(os.path.join(self.folder_path, 'ii.itp')),
            luigi.LocalTarget(os.path.join(self.folder_path, 'out_dg.dat')),
            luigi.LocalTarget(os.path.join(self.folder_path, 'averageA.gro'))
            ]
        for s in self.states:
            targets.append([
                luigi.LocalTarget(os.path.join(self.folder_path, 'dump%s.gro'%s)),
                luigi.LocalTarget(os.path.join(self.folder_path, 'all_eq%s_fit.xtc'%s))
                ])
        for i in range(self.study_settings['n_repeats']):
            for m in range(self.study_settings['n_sampling_sims']):
                targets.append(luigi.LocalTarget(os.path.join(self.folder_path,
                                              "topolTI_ions%d_%d.top"%(i,m))))
        return targets
