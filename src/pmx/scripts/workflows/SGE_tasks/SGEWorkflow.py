#!/usr/bin/env python

import logging
import luigi
from luigi.tools.deps_tree import print_tree
from pmx.scripts.workflows.Workflow import Workflow
from pmx.scripts.workflows.utils import confirm_defNO


# Wrapper task that executes requested dependencies
class SGE_wrapper(luigi.task.WrapperTask):

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="", description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        description="Explicit job name given via qsub.")

    my_deps=[]
    def work(self):
        pass

    def set_deps(self, deps):
        self.my_deps=deps

    def requires(self):
        return( self.my_deps )

# ==============================================================================
#                             Workflow Class
# ==============================================================================
class SGE_Workflow(Workflow):

    def __init__(self, pe='openmp_fast', rem_sched=False,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.states={"A":"l0", "B":"l1"} #states and suffixes of mdp files
        self.TIstates={"A":"l0", "C":"l1"} #states and suffixes of mdp files

        self.tasks=[]
        self.n_workers=1
        self.pe=pe
        self.rem_sched=rem_sched
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
                             'mdrun_opts':self.mdrun_opts,
                             'b':self.b}

        #set luigi logger level to mimimize clutter
        #logger = logging.getLogger('luigi-interface')
        #logger.setLevel(logging.WARNING)

        #sanity checks
        self.check_sanity()
        self.check_inputs()

    def run_everything(self):
        """Runs the whole workflow.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        """
        #Run the tasks
        test=SGE_wrapper()
        test.set_deps(self.tasks)

        print(print_tree(test))
        #exit(1)


        if(self.rem_sched): #use remote scheduler
            mylogger = logging.getLogger(self.__class__.__name__)
            mylogger.warning("Luigi's central scheduler uses an UNSECURED "
                  "protocol (http), so make sure the port it communicates "
                  "on is otherwize secured (eg. IP filtering).")

            if(not confirm_defNO("Do you want to continue anyway?")):
                exit(0)
            else:
                luigi.build([test],
                    local_scheduler=False, workers=self.n_workers)
        else: #use local scheduler
            luigi.build([test],
                    local_scheduler=True, workers=self.n_workers)


