#!/usr/bin/env python

import logging
from pmx.scripts.workflows.utils import NoMissingModuleFilter
if __name__ == '__main__': #mute extraneous missing module warnings from luigi
    logger = logging.getLogger('luigi-interface')
    logger.addFilter(NoMissingModuleFilter())
import luigi
import os
from luigi.tools.deps_tree import print_tree
from pmx.scripts.workflows.utils import parse_options
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.TI import Task_WL_TI_simArray
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.Workflow_aligned_in_Protein import SGE_Workflow_aligned_in_Protein, my_WorkerSchedulerFactory, SGE_test

# ==============================================================================
#                             Workflow Class
# ==============================================================================
class SGE_Workflow_in_Water(SGE_Workflow_aligned_in_Protein):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #overwrite default values
        self.TIstates=self.states #states and suffixes of mdp files

    def gen_folder_name(self,host,ligand):
        return(self.basepath+'/'+host+'/lig_'+ligand)

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
                             'mdrun_opts':self.mdrun_opts,
                             'b':self.b}
        self.tasks=[]

        #TI with no looping over p
        p="water"
        for l in self.ligands:
            folder_path = self.gen_folder_name(p,l)
            for sTI in self.TIstates:
                for i in range(self.n_repeats):
                    for m in range(self.n_sampling_sims):
                        self.tasks.append(Task_WL_TI_simArray(
                            l = l, i = i, m = m, sTI = sTI,
                            study_settings = self.study_settings,
                            folder_path = folder_path,
                            parallel_env='openmp_fast'))

        test=SGE_test()
        test.set_deps(self.tasks)

        print(print_tree(test))
        #exit(1)

        #run SGE_test on login node to bypass scheduler
        n_workers=len(self.hosts)*len(self.ligands)*len(self.states)*\
                    self.n_repeats*self.n_sampling_sims
        luigi.build([test],
                    worker_scheduler_factory=my_WorkerSchedulerFactory(),
                    local_scheduler=True, workers=n_workers)
        #luigi.build([test], workers=2)



# ==============================================================================
#                               MAIN
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

    w=SGE_Workflow_in_Water(toppath, mdppath, ["water"], ["lig"],
                         basepath=basepath,
                         #b=args.b,
                         b=0,
                         #mdrun="mdrun_threads_AVX2_256",
                         mdrun="gmx mdrun",
                         #mdrun_opts="-pin on -nsteps 1000"
                         mdrun_opts="-pin on"
                         )

    w.run_everything()

    print("Complete.\n")

if __name__ == '__main__':
    args = parse_options()
    main(args)