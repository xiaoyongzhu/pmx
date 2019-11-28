import glob
import luigi
import os
import shutil as sh
import datetime
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster

class SGE_Sim(SGETunedJobTask):

    #Parameters:
    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    #values to be overwritten by subclasses:
    base_path = "./"
    sim_path = "./"
    mdp = "input.mdp"
    top = "topol.top"
    struct = "input.gro"
    posre = None
    mdrun = "gmx mdrun"
    mdrun_opts = ""

    def output(self):
        """
        Returns the target output for this task.
        In this case, a successful execution of this task will create a file on the local filesystem.
        :return: the target output for this task.
        :rtype: object (:py:class:`luigi.target.Target`)
        """
        return luigi.LocalTarget(os.path.join(self.sim_path, 'confout.gro'))

    #Override requires() in Workflow specific implementations

    def work(self):
        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)

        #restraint
        r=""
        if self.posre:
            r="-r "+self.posre

        #make tpr (if it wasn't made in past crashed run)
        if(not os.path.isfile("tpr.tpr")):
            os.system("gmx grompp -p %s -c %s %s "
                      "-o tpr.tpr -f %s -v -maxwarn 3 "
                      "> prep.log 2>&1"%(
                          self.top, self.struct, r, self.mdp)
                      )

        #don't run mdrun on failure
        if(not os.path.isfile("tpr.tpr")):
            os.chdir(self.base_path)
            raise Exception("Failed generating file: " +
                            os.path.join(self.sim_path,"tpr.tpr") )

        #run sim
        if(os.path.isfile('state.cpt')):
            #checkpoint exists, resume from it
            date_time = datetime.datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
            sh.copy('state.cpt','state_{}_start.cpt'.format(date_time))
            os.system("{mdrun} -s tpr.tpr -ntomp {n_cpu} -cpi state.cpt "
                      "{mdrun_opts} > mdrun.log 2>&1".format(
                          mdrun=self.mdrun, n_cpu=self.n_cpu,
                          mdrun_opts=self.mdrun_opts))
        else: #first time
            os.system("{mdrun} -s tpr.tpr -ntomp {n_cpu} "
                      "{mdrun_opts} > mdrun.log 2>&1".format(
                          mdrun=self.mdrun, n_cpu=self.n_cpu,
                          mdrun_opts=self.mdrun_opts))

        #clean overwritten files
        cleanList = glob.glob(self.sim_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)

        #Return to basepath
        os.chdir(self.base_path)
