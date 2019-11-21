#!/usr/bin/env python

import logging
from pmx.scripts.workflows.utils import NoMissingModuleFilter
if __name__ == '__main__': #mute extraneous missing module warnings from luigi
    logger = logging.getLogger('luigi-interface')
    logger.addFilter(NoMissingModuleFilter())
import os
from pmx.scripts.workflows.utils import parse_options
from pmx.scripts.workflows.SGE_tasks.SGEWorkflow import SGE_Workflow
from pmx.scripts.workflows.SGE_tasks.absFE.summary import Task_summary_aligned

# ==============================================================================
#                             Workflow Class
# ==============================================================================
class SGE_Workflow_aligned_complete(SGE_Workflow):

    def run_everything(self):
        """Runs the whole workflow.

        Parameters
        ----------
        None.

        Returns
        -------
        None.
        """

        self.tasks.append(Task_summary_aligned(
            hosts = self.hosts, ligands = self.ligands,
            study_settings = self.study_settings,
            parallel_env=self.pe))

        self.n_workers=len(self.hosts)*len(self.ligands)*len(self.states)*\
                    self.n_repeats*self.n_sampling_sims

        super().run_everything() #creates the scheduler and runs the workers




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

    w=SGE_Workflow_aligned_complete(
            toppath=toppath, mdppath=mdppath,
            hosts=["BRD1"],
            ligands=["lig"],
            basepath=basepath,
            #b=args.b,
            b=0,
            mdrun=args.mdrun,
            mdrun_opts=args.mdrun_opts,
            pe=args.pe,
            rem_sched=args.rem_sched
            )

    w.run_everything()

    print("Complete.\n")

if __name__ == '__main__':
    args = parse_options(SGE=True)
    main(args)