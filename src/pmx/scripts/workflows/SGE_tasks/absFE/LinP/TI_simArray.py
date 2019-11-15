import glob
import luigi
import os
import subprocess
import pmx.scripts.workflows.SGE_tasks.SGETunedRunner as sge_runner
from luigi.contrib.sge import logger
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.utils import read_from_mdp



def _parse_qsub_job_array_id(qsub_out):
    """Parse job id from qsub output string.

    Assume format:

        "Your job-array <job_id>.<min>-<max>:<step> ("<job_name>") has been submitted"

    """
    return int(qsub_out.split()[2].split('.')[0])

def _build_qsub_array_command(cmd, job_name, pe, n_cpu, ndHdl_left):
    """Submit array job shell command to SGE queue via `qsub`"""
    qsub_template = """echo {cmd} | qsub -t 1:{ndHdl_left} -V -r y -pe {pe} {n_cpu} -N {job_name} -cwd """
    return qsub_template.format(
        cmd=cmd, job_name=job_name,
        pe=pe, n_cpu=n_cpu, ndHdl_left=ndHdl_left)


class SGE_TI_simArray(SGETunedJobTask):

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

    not_finished = luigi.ListParameter(description='List of not finished '
                      'frame numbers', significant=True,
                      default=[])

    stage="morphes"
    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/ti_{0}.mdp".format(
                self.study_settings['states'][self.sTI])
        self.top = self.folder_path+"/topolTI_ions{3}_{4}.top".format(
            self.p, self.l, self.s, self.i, self.m)

        self.unfinished=[]

    def output(self):
        nframes = len(glob.glob1(self.sim_path,"frame*.gro"))

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'dHdl%d.xvg'%nf)) )
        return targets

    def find_unfinished_dHdl(self):
        expected_last_time = read_from_mdp(self.mdp)
        nframes = len(glob.glob1(self.sim_path,"frame*.gro"))

        self.unfinished=[]
        for nf in range(nframes):
            fname=os.path.join(self.sim_path, 'dHdl%d.xvg'%nf)
            if(os.path.isfile(fname)):
                last_line = subprocess.check_output(
                    ['tail', '-1', fname]).decode('utf-8')
                last_time = float(last_line.split()[0])
                if(last_time == expected_last_time):
                    next;

            #frame not ready
            self.unfinished.append(nf)



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
        job_str = '"source ~/.luigi_profile; '+job_str+'"'

        #find list of unfinished dHdl ids, this should get pickled
        self.find_unfinished_dHdl()

        # Build qsub submit command
        submit_cmd = _build_qsub_array_command(job_str, self.task_family,
                                               self.parallel_env, self.n_cpu,
                                               len(self.unfinished))
        logger.debug('qsub command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True).decode('utf-8')
        logger.debug("Submitted job to qsub with response:\n" + output)
        self.job_id = _parse_qsub_job_array_id(output)
        #logger.debug("Submitted job to qsub with response:\n" + output)

        self._track_job()

        # Now delete the temporaries, if they're there.
        if (self.tmp_dir and os.path.exists(self.tmp_dir) and not self.dont_remove_tmp_dir):
            logger.info('Removing temporary directory %s' % self.tmp_dir)
            subprocess.call(["rm", "-rf", self.tmp_dir])



    def work(self):
        os.makedirs(self.sim_path, exist_ok=True)
        os.chdir(self.sim_path)

        #find which task this is
        SGE_TASK_ID = os.environ['SGE_TASK_ID']

        #translate it into a non-finished frame ID
        dHdL_id=self.unfinished[SGE_TASK_ID]

        #find temp working dir for this job
        TMPDIR = os.environ['TMPDIR']

        #make tpr
        os.system("gmx grompp -p %s -c %s "
                  "-o %s/ti.tpr -f %s -v -maxwarn 3 "%(
                      TMPDIR, self.top, "frame%d.gro"%dHdL_id, self.mdp)
                  )

        #run sim
        os.system(self.mdrun+" -s {D}/ti.tpr -dhdl {D}/dgdl.xvg -cpo "
                  "{D}/state.cpt -e {D}/ener.edr -g {D}/md.log -o "
                  "{D}/traj.trr -x {D}/traj.xtc -c {D}/confout.gro "
                  "-ntomp {n_cpu} {opts}".format(
                      D=TMPDIR, n_cpu=self.n_cpu, opts=self.mdrun_opts)
                  )

        #copy dHdl file back
        os.system("rsync {}/dgdl.xvg {}/dHdl{}.xvg".format(
                      TMPDIR, self.sim_path, dHdL_id)
                  )

        #Return to basepath
        os.chdir(self.base_path)
