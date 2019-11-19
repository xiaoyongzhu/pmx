import errno
import glob
import logging
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

def _build_qsub_array_command(cmd, job_name, pe, n_cpu, work_dir, ndHdl_left):
    """Submit array job shell command to SGE queue via `qsub`"""
    qsub_template = """echo {cmd} | qsub -t 1:{ndHdl_left} -V -r y -pe {pe} {n_cpu} -N {job_name} -wd {work_dir}"""
    return qsub_template.format(
        cmd=cmd, job_name=job_name, work_dir=work_dir,
        pe=pe, n_cpu=n_cpu, ndHdl_left=ndHdl_left)


class Task_TI_simArray(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    sTI = luigi.Parameter(description='Coupling state')

    folder_path = luigi.Parameter(significant=False,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
                 description='Dict of study stettings '
                      'used to propagate settings to dependencies')

    stage="morphes"
    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False)

    job_name_format = luigi.Parameter(
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/ti_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])
        self.top = self.folder_path+"/topolTI_ions{3}_{4}.top".format(
            self.p, self.l, self.sTI, self.i, self.m)
        self.mdrun = self.study_settings['mdrun']
        self.mdrun_opts = self.study_settings['mdrun_opts']

        self.unfinished=[]
        #find list of unfinished dHdl ids, this should get pickled
        self.unfinished=self.find_unfinished_dHdl()

    def output(self):
        nframes = len(glob.glob1(self.sim_path,"frame*.gro"))

        targets=[]
        for nf in range(nframes):
            targets.append(luigi.LocalTarget(
                os.path.join(self.sim_path, 'dHdl%d.xvg'%nf)) )
        return targets

    def complete(self):
        """
        Check if all dHdl files exist and are of correct length.
        """
        outputs = luigi.task.flatten(self.output())
        all_exist=all(map(lambda output: output.exists(), outputs))

        return (all_exist and not self.find_unfinished_dHdl())


    def find_unfinished_dHdl(self):
        expected_end_time, dtframe = read_from_mdp(self.mdp)
        nframes = len(glob.glob1(self.sim_path,"frame*.gro"))

        unf=[]
        for nf in range(nframes):
            fname=os.path.join(self.sim_path, 'dHdl%d.xvg'%nf)
            if(os.path.isfile(fname)):
                last_line = subprocess.check_output(
                    ['tail', '-1', fname]).decode('utf-8')
                last_time = float(last_line.split()[0])
                if(last_time == expected_end_time):
                    continue;

            #frame not ready
            unf.append(nf)
        return(unf)



    def _run_job(self):

        # Build a qsub argument that will run sge_runner.py on the directory we've specified
        runner_path = sge_runner.__file__
        if runner_path.endswith("pyc"):
            runner_path = runner_path[:-3] + "py"
        job_str = 'python {0} "{1}" "{2}"'.format(
            runner_path, self.tmp_dir, os.getcwd())  # enclose tmp_dir in quotes to protect from special escape chars
        if self.no_tarball:
            job_str += ' --no-tarball'

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

        #tell runner that this is an array job
        job_str += ' --arrayjob'

        #force loading of conda and luigi by sourcing a custom profile
        job_str = '"source ~/.luigi_profile; '+job_str+'"'

        self.errfile=""; #no errorfile. mdrun dumps too much into stderr

        # Build qsub submit command
        submit_cmd = _build_qsub_array_command(job_str, self.job_name,
                                               self.parallel_env, self.n_cpu,
                                               self.sim_path,
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
        SGE_TASK_ID = int(os.environ['SGE_TASK_ID'])

        #translate it into a non-finished frame ID
        print("unfinished:",self.unfinished)
        dHdL_id=self.unfinished[SGE_TASK_ID-1] #SGE starts counting from 1

        #find temp working dir for this job
        TMPDIR = os.environ['TMPDIR']

        #make tpr
        os.system("gmx grompp -p {top} -c frame{_id}.gro "
                  "-o {D}/ti.tpr -po {D}/mdout.mdp -f {mdp} "
                  "-v -maxwarn 3 ".format(D=TMPDIR, top=self.top,
                                          mdp=self.mdp, _id=dHdL_id) )

        #run sim
        os.system(self.mdrun+" -s {D}/ti.tpr -dhdl {D}/dgdl.xvg -cpo "
                  "{D}/state.cpt -e {D}/ener.edr -g {D}/md.log -o "
                  "{D}/traj.trr -x {D}/traj.xtc -c {D}/confout.gro "
                  "-ntomp {n_cpu} {opts}".format(
                      D=TMPDIR, n_cpu=self.n_cpu, opts=self.mdrun_opts) )

        #copy dHdl file back
        os.system("rsync {}/dgdl.xvg {}/dHdl{}.xvg".format(
                      TMPDIR, self.sim_path, dHdL_id) )

        #Return to basepath
        os.chdir(self.base_path)
