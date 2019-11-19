import errno
import logging
import luigi
import os
import pmx
import random
import subprocess
import sys
import time
import pmx.scripts.workflows.SGE_tasks.SGETunedRunner as sge_runner
from luigi.contrib.hadoop import create_packages_archive
from luigi.contrib.sge import SGEJobTask, logger, _parse_qstat_state, _build_qsub_command, _parse_qsub_job_id
try:
    import cPickle as pickle
except ImportError:
    import pickle


class SGETunedJobTask(SGEJobTask):

#TODO: Override _track_job() so that:
#           - support batching for TI
#           - need resume support in case job runs out of time

    #change default parallel environment
    parallel_env = luigi.Parameter(default='openmp_fast', significant=False)
    #poll time in seconds.
    #Large value is better to avoid overloading the login node.
    #Needs to be less than time between queue updates.
    poll_time = luigi.IntParameter(
        significant=False, default=30,
        description="specify the wait time to poll qstat for the job status")

    #temp files
    shared_tmp_dir = luigi.Parameter(default=os.path.join(os.getenv("HOME"), 'temp'), significant=False)
    # dont_remove_tmp_dir = luigi.BoolParameter(
    #     significant=False,
    #     default=True,
    #     description="don't delete the temporary directory used (for debugging)")

    #Don't archive luigi or pmx. Jobs load them through conda
    # no_tarball = luigi.BoolParameter(
    #     significant=False,
    #     default=True,
    #     description="don't tarball (and extract) the luigi project files")

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}", description="A string that can be "
        "formatted with class variables to name the job with qsub.")
    job_name = luigi.Parameter(
        significant=False, default="",
        description="Explicit job name given via qsub.")

    disable_window=3600*24*7 # 7 days
    retry_count=0 #no retries within disable_window seconds of previous failure

    extra_packages=[] #extra packages to be tarballed. Overloaded by subclasses.

    def _init_local(self):

        # Set up temp folder in shared directory (trim to max filename length)
        base_tmp_dir = self.shared_tmp_dir
        random_id = '%016x' % random.getrandbits(64)
        folder_name = self.task_id + '-' + random_id
        self.tmp_dir = os.path.join(base_tmp_dir, folder_name)
        max_filename_length = os.fstatvfs(0).f_namemax
        self.tmp_dir = self.tmp_dir[:max_filename_length]
        logger.info("Tmp dir: %s", self.tmp_dir)
        os.makedirs(self.tmp_dir)

        # Dump the code to be run into a pickle file
        logging.debug("Dumping pickled class")
        self._dump(self.tmp_dir)

        if not self.no_tarball:
            # Make sure that all the class's dependencies are tarred and available
            # This is not necessary if luigi is importable from the cluster node
            logging.debug("Tarballing dependencies")
            # Grab luigi, the whole of pmx, and the module containing the code to be run
            packages = [luigi] + [pmx] + self.extra_packages +\
                [__import__(self.__module__, None, None, 'dummy')]
            create_packages_archive(packages, os.path.join(self.tmp_dir, "packages.tar"))

    def _run_job(self):

        # Build a qsub argument that will run sge_runner.py on the directory we've specified
        runner_path = sge_runner.__file__
        if runner_path.endswith("pyc"):
            runner_path = runner_path[:-3] + "py"
        job_str = 'python {0} "{1}" "{2}"'.format(
            runner_path, self.tmp_dir, os.getcwd())  # enclose tmp_dir in quotes to protect from special escape chars
        if self.no_tarball:
            job_str += ' "--no-tarball"'

            #force loading of dependencies by sourcing a custom profile
            if(os.path.isfile("~/.luigi_profile")):
                job_str = '"source ~/.luigi_profile; '+job_str+'"'
            else:
                logging.error("Tarballing of dependencies is disabled and "
                              "~/.luigi_profile does not exist. "
                              "Will not be able to load all workflow "
                              "dependencies without it. Please create it and "
                              "within activate a conda environment containing "
                              "at least python>3.6, "
                              "pmx, luigi, MDanalysis, matplotlib, and numpy.")
                raise FileNotFoundError(errno.ENOENT,
                              os.strerror(errno.ENOENT), "~/.luigi_profile")

        # Build qsub submit command
        self.outfile = os.path.join(self.tmp_dir, 'job.out')
        self.errfile = os.path.join(self.tmp_dir, 'job.err')
        submit_cmd = _build_qsub_command(job_str, self.job_name, self.outfile,
                                         self.errfile, self.parallel_env, self.n_cpu)
        logger.debug('qsub command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True).decode('utf-8')
        logger.debug("Submitted job to qsub with response:\n" + output)
        self.job_id = _parse_qsub_job_id(output)
        #logger.debug("Submitted job to qsub with response:\n" + output)

        self._track_job()

        # Now delete the temporaries, if they're there.
        if (self.tmp_dir and os.path.exists(self.tmp_dir) and not self.dont_remove_tmp_dir):
            logger.info('Removing temporary directory %s' % self.tmp_dir)
            subprocess.call(["rm", "-rf", self.tmp_dir])

    def _track_job(self):
        while True:
            # Sleep for a little bit
            time.sleep(self.poll_time)

            # See what the job's up to
            # ASSUMPTION
            qstat_out = subprocess.check_output(['qstat']).decode('utf-8')
            sge_status = _parse_qstat_state(qstat_out, self.job_id)
            if sge_status == 'r':
                logger.info('Job is running...')
            elif sge_status == 'qw':
                logger.info('Job is pending...')
            elif sge_status == 't':
                logger.info('Job is transferring...')
            elif 'E' in sge_status:
                logger.error('Job has FAILED:\n' + '\n'.join(self._fetch_task_failures()))
                break
            elif sge_status == 'u':
                # Then the job could either be failed or done.
                errors = self._fetch_task_failures()
                if not errors:
                    logger.info('Job is done')
                else:
                    logger.error('Job has FAILED:\n' + '\n'.join(errors))
                break
            else:
                logger.info('Job status is UNKNOWN!')
                logger.info('Status is : %s' % sge_status)
                raise Exception("job status isn't one of ['r', 'qw', 'E*', 't', 'u']: %s" % sge_status)

    def _dump(self, out_dir=''):
        """Dump instance to file."""
        with self.no_unpicklable_properties():
            self.job_file = os.path.join(out_dir, 'job-instance.pickle')
            if self.__module__ == '__main__':
                d = pickle.dumps(self)
                module_name = os.path.basename(sys.argv[0]).rsplit('.', 1)[0]
                d = d.replace(b'c__main__', b"c" + module_name.encode('utf-8'))
                open(self.job_file, "wb").write(d)
            else:
                pickle.dump(self, open(self.job_file, "wb"))

