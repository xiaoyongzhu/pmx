import os
import subprocess
import time
from luigi.contrib import sge_runner
from luigi.contrib.sge import SGEJobTask, logger, _parse_qstat_state, _build_qsub_command, _parse_qsub_job_id

class SGETunedJobTask(SGEJobTask):

#TODO: Override _track_job() so that:
#           - It looks for the specific job id not the whole list of jobs
#           - sge_status == 't' is not nessesarily done
#           - check if job has left the queue
#           - support batching for TI
# see https://luigi.readthedocs.io/en/stable/_modules/luigi/contrib/sge.html#SGEJobTask._track_job

    def _run_job(self):

        # Build a qsub argument that will run sge_runner.py on the directory we've specified
        runner_path = sge_runner.__file__
        if runner_path.endswith("pyc"):
            runner_path = runner_path[:-3] + "py"
        job_str = 'python {0} "{1}" "{2}"'.format(
            runner_path, self.tmp_dir, os.getcwd())  # enclose tmp_dir in quotes to protect from special escape chars
        if self.no_tarball:
            job_str += ' "--no-tarball"'

        #force loading of conda and luigi by sourcing a custom profile
        job_str = "source ~/.luigi_profile; "+job_str

        # Build qsub submit command
        self.outfile = os.path.join(self.tmp_dir, 'job.out')
        self.errfile = os.path.join(self.tmp_dir, 'job.err')
        submit_cmd = _build_qsub_command(job_str, self.task_family, self.outfile,
                                         self.errfile, self.parallel_env, self.n_cpu)
        logger.debug('qsub command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True)
        self.job_id = _parse_qsub_job_id(output)
        logger.debug("Submitted job to qsub with response:\n" + output)

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
            qstat_out = subprocess.check_output(['qstat'])
            sge_status = _parse_qstat_state(qstat_out, self.job_id)
            if sge_status == 'r':
                logger.info('Job is running...')
            elif sge_status == 'qw':
                logger.info('Job is pending...')
            elif 'E' in sge_status:
                logger.error('Job has FAILED:\n' + '\n'.join(self._fetch_task_failures()))
                break
            elif sge_status == 't' or sge_status == 'u':
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