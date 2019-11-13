#!/usr/bin/env python

import glob
import logging
import luigi
import os
import shutil as sh
from luigi.contrib.sge import LocalSGEJobTask
from pmx.scripts.workflows.Workflow import check_file_ready
from pmx.scripts.workflows.Workflow_alligned_in_protein import Workflow_alligned_inProtein, parse_options
from pmx.scripts.workflows.SGE_tasks.Sim import SGE_Sim

# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Gather_Inputs_PL_folder(LocalSGEJobTask): # will execute on the login node
    folder_path = luigi.Parameter(significant=False,
                      description='Path to the protein+ligand folder to set up')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

    study_settings = luigi.DictParameter(description='Dict of study stettings '
                      'used to propagate settings to dependencies')

    def work(self):

        #make folder
        os.makedirs(self.folder_path, exist_ok=True)
        os.chdir(self.folder_path)

        #topology
        sh.copy(self.study_settings['top_path']+"/topol_abs_prot_norestr_amber.top",
                self.folder_path+"/topol.top")
        sh.copy(self.study_settings['top_path']+"/ligand/"+self.l+"/lig.itp",self.folder_path+"/lig.itp")
        sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot.itp",self.folder_path+"/prot.itp")

        #initial coordinates where protein and ligand are bound
        sh.copy(self.study_settings['top_path']+"/proteins/"+self.p+"/prot_lig.pdb",
                self.folder_path+"/init.pdb")

        #generate temporary index file
        os.system("echo 'q\n' | gmx make_ndx -f init.pdb "
                  "-o index.ndx > setup.log 2>&1")
        check_file_ready("index.ndx")

        #generate restraints for equillibration
        #TODO: rewrite this to use the pmx Topology class
        os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                  "-fc 9000 9000 9000 -o prot_posre.itp "
                  "-n index.ndx >> setup.log 2>&1")
        check_file_ready("prot_posre.itp")
        os.system("echo 'Protein\n' | gmx genrestr -f init.pdb "
                  "-fc 500 500 500 -o prot_posre_soft.itp "
                  "-n index.ndx >> setup.log 2>&1")
        check_file_ready("prot_posre_soft.itp")
        os.system("echo 'MOL\n' | gmx editconf -f init.pdb "
                  "-o lig.pdb -n index.ndx >> setup.log 2>&1")
        check_file_ready("lig.pdb")
        os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                  "-fc 9000 9000 9000 "
                  "-o lig_posre.itp >> setup.log  2>&1")
        check_file_ready("lig_posre.itp")
        os.system("echo 'MOL\n' | gmx genrestr -f lig.pdb "
                  "-fc 500 500 500 "
                  "-o lig_posre_soft.itp >> setup.log 2>&1")
        check_file_ready("lig_posre_soft.itp")

        #clean overwritten files
        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                raise OSError("Error while deleting file: "+filePath)

        #Return to basepath
        os.chdir(self.study_settings['base_path'])

    def output(self):
        files=["topol.top", "lig.itp", "prot.itp", "init.pdb", "index.ndx",
               "prot_posre.itp", "prot_posre_soft.itp",
               "lig_posre.itp", "lig_posre_soft.itp"]
        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


class Prep_PL_folder(LocalSGEJobTask): # will execute on the login node
    folder_path = luigi.Parameter(significant=False,
                         description='Path to the protein+ligand folder to set up')
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')

    study_settings = luigi.DictParameter(description='Dict of study stettings '
                      'used to propagate settings to dependencies')

    def requires(self):
        return( Gather_Inputs_PL_folder(folder_path=self.folder_path,
                                        p=self.p, l=self.l,
                                        study_settings=self.study_settings) )

    def work(self):
        os.chdir(self.folder_path)

        os.system("gmx editconf -f init.pdb -o box.pdb -bt %s -d %f "\
                  "> prep.log 2>&1"%(self.study_settings['bt'], self.study_settings['d']))
        check_file_ready("box.pdb")
        sh.copy("topol.top","topol_solvated.top")
        os.system("gmx solvate -scale 1.0 -cp box.pdb -o water.pdb "\
                  "-cs spc216.gro -p topol_solvated.top >> prep.log 2>&1")
        check_file_ready("water.pdb")
        os.system("gmx grompp -p topol_solvated.top -c water.pdb -o tpr.tpr "\
                  "-f %s/protein/init.mdp -v -maxwarn 2 "\
                  ">> prep.log 2>&1"%self.study_settings['mdp_path'])
        check_file_ready("tpr.tpr")

        #generate ions for each
        #independent repeat (multiple for confidence estimate)
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                top_ions="topol_ions%d_%d.top"%(i,m)
                pdb_ions="ions%d_%d.pdb"%(i,m)
                if(os.path.isfile(pdb_ions)): #skip if it already exists
                    continue
                sh.copy("topol_solvated.top", top_ions)
                os.system("echo 'SOL' | gmx genion -s tpr.tpr "
                          "-p %s -conc %f "
                          "-neutral -nname Cl -pname Na "
                          "-o %s >> genion.log 2>&1" %(
                              top_ions, self.study_settings['salt_conc'], pdb_ions) )
                check_file_ready(pdb_ions)

        cleanList = glob.glob(self.folder_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)

        #Return to basepath
        os.chdir(self.study_settings['base_path'])

    def output(self):
        files=[]
        #independent repeat (multiple for confidence estimate)
        for i in range(self.study_settings['n_repeats']):
            #sampling simulations in each repeat
            for m in range(self.study_settings['n_sampling_sims']):
                files.append("topol_ions%d_%d.top"%(i,m))
                files.append("ions%d_%d.pdb"%(i,m))

        return [luigi.LocalTarget(os.path.join(self.folder_path, f)) for f in files]


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

# ==============================================================================
#                             Workflow Class
# ==============================================================================
class SGE_Workflow_alligned_inProtein(Workflow_alligned_inProtein):

    def run_everything(self):
        """Runs the whole workflow.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        """

        #sanity checks
        self.check_sanity()
        self.check_inputs()

        logger = logging.getLogger('luigi-interface')
        logger.setLevel(logging.WARNING)

        self.study_settings={'base_path':self.basepath,
                             'top_path':self.toppath,
                             'mdp_path':self.mdppath,
                             'd':self.d, 'bt':self.bt,
                             'salt_conc':self.salt_conc,
                             'n_repeats':self.n_repeats,
                             'n_sampling_sims':self.n_sampling_sims,
                             'states':self.states,
                             'TIstates':self.TIstates,
                             'mdrun':self.mdrun,
                             'mdrun_opts':self.mdrun_opts }
        self.tasks=[]

        #run NPT to sample starting frames for TI
        for p in self.hosts:
            for l in self.ligands:
                folder_path = self.gen_folder_name(p,l)
                for s in self.states:
                    for i in range(self.n_repeats):
                        for m in range(self.n_sampling_sims):
                            self.tasks.append(Sim_PL_NPT(
                                p = p, l = l, i = i, m = m, s = s,
                                study_settings = self.study_settings,
                                folder_path = folder_path,
                                parallel_env='openmp_fast'))

        #genergate Boresh-style protein-ligand restraints

        #align vaccum ligand onto apo protein structures

        #TI

        #Run the tasks
        class SGE_test(LocalSGEJobTask): # will execute on the login node

            my_deps=[]
            def work(self):
                pass

            def set_deps(self, deps):
                self.my_deps=deps

            def requires(self):
                return( self.my_deps )

        test=SGE_test()
        test.set_deps([self.tasks[0]])

        #run SGE_test on login node to bypass scheduler
        luigi.build([test], local_scheduler=True, retry_count=0, workers=16)
        #luigi.build([test], workers=2)



# ==============================================================================
#                               FUNCTIONS
# ==============================================================================

def main(args):
    """Run the main script.

    Parameters
    ----------
    args : argparse.Namespace
        The command line arguments
    """
    toppath=os.path.abspath(args.toppath)
    mdppath=os.path.abspath(args.mdppath)
    basepath=os.path.abspath(args.basepath)

    w=SGE_Workflow_alligned_inProtein(toppath, mdppath, ["BRD1"], ["lig"],
                         basepath=basepath,
                         #mdrun="mdrun_threads_AVX2_256",
                         mdrun="gmx mdrun",
                         mdrun_opts="-pin on -nsteps 1000")

    w.run_everything()

    print("Complete.\n")

if __name__ == '__main__':
    args = parse_options()
    main(args)