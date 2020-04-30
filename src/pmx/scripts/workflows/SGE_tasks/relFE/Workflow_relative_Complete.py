#!/usr/bin/env python

import logging
from pmx.scripts.workflows.utils import NoMissingModuleFilter
if __name__ == '__main__': #mute extraneous missing module warnings from luigi
    logger = logging.getLogger('luigi-interface')
    logger.addFilter(NoMissingModuleFilter())
import os
from pmx.scripts.workflows.utils import parse_options
from pmx.scripts.workflows.SGE_tasks.SGEWorkflow import SGE_Workflow
from pmx.scripts.workflows.SGE_tasks.relFE.summary import Task_summary_relative

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
        summary=Task_summary_relative(
            #hosts = self.hosts, ligands = self.ligands,
            study_settings = self.study_settings,
            parallel_env=self.pe)
        self.tasks.append(summary)

        self.n_workers=2*len(summary.hosts)*len(summary.ligands)*len(self.states)*\
                    summary.n_repeats*summary.n_sampling_sims

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
            basepath=basepath,
            b=args.b,
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
