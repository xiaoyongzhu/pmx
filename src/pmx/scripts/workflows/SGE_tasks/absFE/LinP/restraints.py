#!/usr/bin/env python

import luigi
import MDAnalysis as md
import os
import sys
import glob
from io import StringIO
from luigi.parameter import ParameterVisibility
from pmx import ndx
from pmx.forcefield import Topology
from pmx.ndx import IndexGroup
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.find_anchors_and_write_ii import find_restraints
from pmx.scripts.workflows.find_anchors_and_write_ii_single_traj import find_restraints_align2crystal
from pmx.scripts.workflows.find_avg import find_avg_struct
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.equil_sims import Sim_PL_NPT
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.alignment import Task_PL_align
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.postHoc_restraining_python3 import main as main_postHock_restr

# ==============================================================================
#                            Helper Functions
# ==============================================================================
def clean_virtual_sites_from_ndx(ndx_fn, mol_sel, filteree, itp):
    my_ndx=ndx.IndexFile(ndx_fn, verbose=False)
    whole_mol_ids = my_ndx[mol_sel].ids
    filteree_ids = my_ndx[filteree].ids
    first_at = whole_mol_ids[0]
    
    top = Topology(itp, is_itp=True, assign_types=False)
    vsites=[] #indeces of virtual sites in ligand
    for l in [top.virtual_sites2, top.virtual_sites3, top.virtual_sites4]:
        for vs in l:
            vsites.append(vs[0].id - 1 + first_at)
    
    new_filteree_ids = []
    for i in filteree_ids:
        if not (i in vsites):
            new_filteree_ids.append(i)
    g = IndexGroup(ids=new_filteree_ids, name=filteree+"_&_!vsites")
    my_ndx.add_group(g)
    my_ndx.write(ndx_fn)

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_gen_restraints(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')

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
        
    #debug output
    debug = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default=False,
        description="show debug output in a log.")

    extra_packages=[md]
    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']

        if(self.restr_scheme=="Aligned"):
            self.states=["C"]
        elif(self.restr_scheme=="Fitted"):
            self.states=["A", "B"]
        elif(self.restr_scheme=="Fixed"):
            raise(Exception("Fixed restraints not yet implemened."))
            self.states = self.study_settings['states']
        else:
            raise(Exception("Unrecognized restraint scheme '%s'"%self.restr_scheme))

    def work(self):
        #generate morphs for A state
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        # srctpr=self.folder_path+"/state{2}/repeat{3}/npt{4}/tpr.tpr"
        # srctraj=self.folder_path+"/state{2}/repeat{3}/npt{4}/traj.trr"

        # # #create prot+MOL index group
        # # os.system("echo \"1|13\nq\n\" | "
        # #           "gmx make_ndx -f ions0_0.pdb "
        # #           "-o index_prot_mol.ndx > /dev/null 2>&1")

        # #make topology for prot+MOL
        # os.system("sed 's/SOL/;SOL/g' topol.top > topol_prot_mol_{i}.top".format(i=self.i))

        # #make topology for ApoP
        # os.system("sed 's/SOL/;SOL/g' {base}/prot_{p}/apoP/topol.top "
        #           "> topol_prot_{i}.top".format(base=self.base_path, p=self.p, i=self.i))

        # for s in self.states:
        #     #make tprs
        #     if(s == "A" or s=="C"):   #align A/C to initial structure
        #         ref="box.pdb"
        #         ref_top="topol_prot_mol_{i}.top".format(i=self.i)
        #         mdp = self.study_settings['mdp_path'] + "/protein/init.mdp"
        #     else:           #align B to average of A
        #         ref="averageA_{i}_prot_only.gro".format(i=self.i)
        #         ref_top="topol_prot_{i}.top".format(i=self.i)
        #         mdp = self.study_settings['mdp_path'] + "/apo_protein/init.mdp"


        #     os.system("gmx grompp -p {ref_top} -c {ref} -f {mdp} "
        #               "-o tpr{s}_{i}.tpr -maxwarn 2 > grompp{s}_{i}.log 2>&1".format(
        #                   ref_top=ref_top, ref=ref, mdp=mdp, s=s, i=self.i) )
        #     check_file_ready("tpr{s}_{i}.tpr".format(s=s, i=self.i))

        #     #collect trjs
        #     # print("\tCollecting trajectories for state%s"%s)

        #     #remove previous log if it exists from a crashed attempt
        #     if(os.path.isfile("trjconv_{i}.log".format(i=self.i))):
        #         os.unlink("trjconv_{i}.log".format(i=self.i))

        #     #independent repeats for error analysis
        #     local_trjs=""

        #     for m in range(self.study_settings['n_sampling_sims']):
        #         if(s=="ApoProt"):
        #             apoP_path=self.base_path+"/prot_{0}/apoP/repeat{3}/npt{4}/"
        #             tpr=apoP_path.format(self.p,self.l,s,self.i,m)+"tpr.tpr"
        #             trj=apoP_path.format(self.p,self.l,s,self.i,m)+"traj.trr"
        #             sel="4 Protein"
        #             ndx=self.base_path+"/prot_{0}/apoP/index.ndx".format(self.p)
        #         if(s=="C"): #explisitly simulated stateC (self.restr_scheme=="Fitted")
        #             aligned_path=self.folder_path+"/state{2}/repeat{3}/{5}{4}/"
        #             tpr=aligned_path.format(self.p,self.l,"A",self.i,m,"npt")+"tpr.tpr"
        #             trj=aligned_path.format(self.p,self.l,s,self.i,m,"morphes")+"aligned.trr"
        #             sel="4 Protein_MOL"
        #             ndx="index_prot_mol.ndx"
        #         else:
        #             tpr=srctpr.format(self.p,self.l,s,self.i,m)
        #             trj=srctraj.format(self.p,self.l,s,self.i,m)
        #             sel="4 Protein_MOL"
        #             ndx="index_prot_mol.ndx"

        #         os.system("echo %s | "
        #                   "gmx trjconv -s %s -f %s "
        #                   "-o eq%s%d_%d.xtc "
        #                   "-sep -ur compact -pbc mol -center "
        #                   "-boxcenter zero -n %s "
        #                   "-b %d >> trjconv.log 2>&1"%(
        #                           sel,tpr,trj, s,self.i,m, ndx,
        #                           self.study_settings['b']) )

        #         check_file_ready("eq%s%d_%d.xtc"%(s,self.i,m))
        #         local_trjs+="eq%s%d_%d.xtc "%(s,self.i,m)

        #     #concatenate trajectories
        #     os.system("gmx trjcat -f {trjs} -o all_eq{s}_{i}.xtc -sort "
        #               "-cat >> trjconv_{i}.log 2>&1".format(
        #                   trjs=local_trjs,s=s,i=self.i) )
        #     check_file_ready("all_eq{s}_{i}.xtc".format(s=s,i=self.i))

        #     #fit to reference structure in tpr files
        #     os.system("echo 4 0 | gmx trjconv -s tpr{s}_{i}.tpr -f all_eq{s}_{i}.xtc "
        #               "-o all_eq{s}_{i}_fit.xtc -fit rot+trans "
        #               ">> trjconv_{i}.log 2>&1".format(s=s,i=self.i) )
        #     check_file_ready("all_eq{s}_{i}_fit.xtc".format(s=s,i=self.i) )

        #     #dump first frame
        #     os.system("echo 0 | gmx trjconv -f all_eq{s}_{i}_fit.xtc "
        #               "-s tpr{s}_{i}.tpr -o dump{s}_{i}.gro -dump 0 "
        #               ">> trjconv_{i}.log 2>&1".format(s=s,i=self.i) )
        #     check_file_ready("dump{s}_{i}.gro".format(s=s,i=self.i))

        #     #find avg structure of A
        #     if(s=="A"):
        #         # print("\tFinding average structure")
        #         find_avg_struct("dumpA_{}.gro".format(self.i), "all_eqA_{i}_fit.xtc".forat(self.i),
        #                         "averageA_{}.gro".format(self.i))
        #         check_file_ready("averageA_{}.gro".format(self.i))
        #         # print("\tExtracting prot. only")

        #         os.system("echo Protein | gmx trjconv -s tprA_{i}.tpr -f averageA_{i}.gro "
        #                   "-o averageA_{i}_prot_only.gro >> trjconv.log 2>&1".format(i=self.i) )
        #         check_file_ready("averageA_{i}_prot_only.gro".format(i=self.i))


        if(self.debug):
            print("debug: restr_scheme={}".format(self.restr_scheme))
        #generate the restraints
        # print("\tGenerating the restraints")
        if(self.restr_scheme=="Aligned"):
            # find_restraints_align2crystal(struct= 'dumpC_{i}.gro'.format(i=self.i),
            #             traj = "all_eqC_{i}_fit.xtc".format(i=self.i),
            #             out="ii_{i}.itp".format(i=self.i),
            #             an_cor_file="out_dg_{i}.dat".format(i=self.i),
            #             plotfile="restraint_coord_distrib_{i}.png".format(i=self.i),
            #             log=False)

            #find correct group indeces
            base_ndx=ndx.IndexFile("index_prot_mol.ndx", verbose=False)
            mol_id = base_ndx.get_group_id("MOL")
            my_ndx_file="index_prot_mol_noH_{i}.ndx".format(i=self.i)
            os.system("echo \"{mol_id} & ! a H*\n\nq\n\" | "
                      "gmx make_ndx -f ions{i}_0.pdb -n index_prot_mol.ndx "
                      "-o {ndx} > noH_make_ndx_{i}.log 2>&1".format(i=self.i, ndx=my_ndx_file, mol_id=mol_id))
            #modify index file to filter out virtual cites from ligand selection
            clean_virtual_sites_from_ndx(my_ndx_file, "MOL", "MOL_&_!H*", "lig.itp")
                      
            check_file_ready(os.path.join(my_ndx_file))
            if(self.debug):
                print("debug: made {}".format(my_ndx_file))
                
                
            aligned_path=self.folder_path+"/state{2}/repeat{3}/{5}{4}/"
            aligned_trjs=""
            for m in range(self.study_settings['n_sampling_sims']):
                aligned_trjs+=aligned_path.format(self.p,self.l,"C",self.i,m,"morphes")+"/frame*.gro"
            
            # os.system("echo -e \"3\n22\n\" | python /home/ykhalak/custom_scripts/pmx/postHoc_restraining_python3.py "
            #           "-f {ap} "
            #           "-n {ndx} -oii ii_{i}.itp -odg out_dg_{i}.dat > gen_restr{i}.log 2>&1".format(
            #               ap=aligned_trjs, ndx=ndx, i=self.i))
                
            if(not os.path.isfile("ii_{i}.itp".format(i=self.i)) or  not os.path.isfile("out_dg_{i}.dat".format(i=self.i))):
                oldstdin = sys.stdin
                oldstdout = sys.stdout
                oldstderr = sys.stderr
                
                my_ndx=ndx.IndexFile(my_ndx_file, verbose=False)
                prot_id = my_ndx.get_group_id("C-alpha")
                mol_id = my_ndx.get_group_id("MOL_&_!H*_&_!vsites")
                sys.stdin = StringIO( "{}\n{}\n".format(prot_id, mol_id) )
                with open("gen_restr{i}.log".format(i=self.i), 'w') as logf:
                    sys.stdout = logf
                    sys.stderr = logf
                    
                    g=glob.glob(aligned_trjs)
                    argv = ["postHoc_restraining_python3.py", "-f", *g, "-n", my_ndx_file,
                                "-oii", "ii_{i}.itp".format(i=self.i),
                                "-odg", "out_dg_{i}.dat".format(i=self.i)]

                    if(self.debug):
                        print("debug: starting postHoc_restraining")                            
                    main_postHock_restr(argv)
                    if(self.debug):
                        print("debug: after postHoc_restraining")

                sys.stdin = oldstdin
                sys.stdout = oldstdout
                sys.stderr = oldstderr

        elif(self.restr_scheme=="Fitted"):
            if(self.debug):
                print("debug: launch find_restraints")
            find_restraints(
                structA="dumpA_{i}.gro".format(i=self.i),
                structB="dumpB_{i}.gro".format(i=self.i),
                ref="averageA_{i}.gro".format(i=self.i),
                trajA="all_eqA_{i}_fit.xtc".format(i=self.i),
                trajB="all_eqB_{i}_fit.xtc".format(i=self.i),
                out="ii_{i}.itp".format(i=self.i),
                an_cor_file="out_dg_{i}.dat".format(i=self.i),
                plotfile="restraint_coord_distrib_{i}.png".format(i=self.i),
                log=False)

        check_file_ready(os.path.join("ii_{i}.itp".format(i=self.i)))

        #create a C state topology that holds ligand in place
        #with the restraint from the ii.itp files

        for m in range(self.study_settings['n_sampling_sims']):
            for s in ["A","C"]:
                nWater,nCl,nNa=(0,0,0)
                top_ions="topol_ions%d_%d.top"%(self.i,m)
                if(s=="C"): # Apo sim can have different number of waters (if the whole protocol wasn't rerun from scratch after including Apo crystal structure)
                    top_Apo=self.folder_path+"/../apoP/topol_ions%d_%d.top"%(self.i,m)
                    with open(top_Apo, 'r') as reftop:
                        for l in reftop:
                            if ("SOL " in l):
                                nWater = int(l.split()[1])
                            if ("Cl " in l):
                                nCl = int(l.split()[1])
                            if ("Na " in l):
                                nNa = int(l.split()[1])
                    
                check_file_ready(top_ions)
                topAC_ions="topolTI_ions%s%d_%d.top"%(s,self.i,m)
                #os.system("cp {} {} > /dev/null 2>&1".format(top_ions,topAC_ions))
                with open(topAC_ions, 'w') as top:
                    with open(top_ions, 'r') as reftop:
                        for l in reftop:
                            if (s=="C" and "SOL " in l):
                                top.write("SOL    {}\n".format(nWater))
                            elif (s=="C" and "Cl " in l):
                                top.write("Cl     {}\n".format(nCl))
                            elif (s=="C" and "Na " in l):
                                top.write("Na     {}\n".format(nNa))
                            else:
                                top.write(l)
                    
                    top.write("\n; Include intermolecular restraints\n")
                    top.write("#include \"ii_{i}.itp\"\n".format(i=self.i))

                #topsrc=self.study_settings['top_path']+"/topol_abs_prot_restr_amber.top"
                #os.system("tail -n 3 {} >> {}".format(topsrc,topAC_ions))
                check_file_ready(topAC_ions)

        #restore base path
        os.chdir(self.base_path)


    def requires(self):
        reqs=[]
        #sampling simulations in each repeat
        for m in range(self.study_settings['n_sampling_sims']):
            #states of equilibrium sims
            if(self.restr_scheme=="Aligned"):
                reqs.append(Task_PL_align(p=self.p, l=self.l, i=self.i, m=m, sTI='C',
                                  study_settings=self.study_settings,
                                  folder_path=self.folder_path,
                                  parallel_env=self.parallel_env,
                                  restr_scheme=self.restr_scheme)
                        )
            elif(self.restr_scheme=="Fitted"):
                for s in self.states:
                    reqs.append(
                        Sim_PL_NPT(p=self.p, l=self.l, i=self.i, m=m, s=s,
                                  study_settings=self.study_settings,
                                  folder_path=self.folder_path,
                                  parallel_env=self.parallel_env)
                        )


        return(reqs)

    def output(self):
        targets=[
            luigi.LocalTarget(os.path.join(self.folder_path, "ii_{i}.itp".format(i=self.i))),
            luigi.LocalTarget(os.path.join(self.folder_path, "out_dg_{i}.dat".format(i=self.i))),
            ]
        if(self.restr_scheme=="Fitted"):
            targets.append([luigi.LocalTarget(os.path.join(self.folder_path, "averageA_{i}.gro".format(i=self.i)))])
        # for s in self.states:
        #     targets.append([
        #         luigi.LocalTarget(os.path.join(self.folder_path, "dump{s}_{i}.gro".format(s=s,i=self.i))),
        #         luigi.LocalTarget(os.path.join(self.folder_path, "all_eq{s}_{i}_fit.xtc".format(s=s,i=self.i)))
        #         ])

        for m in range(self.study_settings['n_sampling_sims']):
            for s in ["A","C"]:
                targets.append(luigi.LocalTarget(os.path.join(self.folder_path,
                                          "topolTI_ions%s%d_%d.top"%(s,self.i,m))))
        return targets
