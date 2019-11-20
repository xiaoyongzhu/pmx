import luigi
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.TI import Task_PL_TI_simArray
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.morphes import Task_WL_gen_morphes


class Task_WL_TI_simArray(Task_PL_TI_simArray):

    #Parameters:
    p = luigi.Parameter(significant=False, default=None,
        description='Protein name') #disables base class' p

    folder_path = luigi.Parameter(significant=False,
                 description='Path to the water+ligand folder to set up')

    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #overwrite relevant variables
        self.mdp = self.study_settings['mdp_path'] +\
            "/water/ti_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])
        self.top = self.folder_path+"/topol_ions{3}_{4}.top".format(
            self.p, self.l, self.sTI, self.i, self.m)
        #self._find_unfinished() is executed later and depends on correct self.mdp

    def requires(self):
        return( Task_WL_gen_morphes(p=self.p, l=self.l,
                          i=self.i, m=self.m, sTI=self.sTI,
                          study_settings=self.study_settings,
                          folder_path=self.folder_path,
                          parallel_env=self.parallel_env) )
