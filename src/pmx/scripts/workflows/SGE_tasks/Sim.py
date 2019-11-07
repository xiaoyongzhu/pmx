import glob
import luigi
import os
from luigi.contrib.sge import SGEJobTask, LocalSGEJobTask

#class SGE_Sim(SGEJobTask):
class SGE_Sim(LocalSGEJobTask): #test locally

    #Parameters:
    study_settings = luigi.DictParameter(description='Dict of study stettings '
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

#TODO: Override _track_job() so that:
#           - It looks for the specific job id not the whole list of jobs
#           - sge_status == 't' is not nessesarily done
#           - check if job has left the queue
#           - support batching for TI
# see https://luigi.readthedocs.io/en/stable/_modules/luigi/contrib/sge.html#SGEJobTask._track_job

    def work(self):
        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)

        #restraint
        r=""
        if self.posre:
            r="-r "+self.posre

        #make tpr
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
        os.system(self.mdrun+" -s tpr.tpr "+\
                  self.mdrun_opts+" > mdrun.log 2>&1")

        #clean overwritten files
        cleanList = glob.glob(self.sim_path+'/#*')
        for filePath in cleanList:
            try:
                os.unlink(filePath)
            except:
                print("Error while deleting file: ", filePath)

        #Return to basepath
        os.chdir(self.base_path)
