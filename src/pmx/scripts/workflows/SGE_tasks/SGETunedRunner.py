# -*- coding: utf-8 -*-
#
# Copyright 2012-2015 Spotify AB
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
The SunGrid Engine runner

The main() function of this module will be executed on the
compute node by the submitted job. It accepts as a single
argument the shared temp folder containing the package archive
and pickled task to run, and carries out these steps:

- extract tarfile of package dependencies and place on the path
- unpickle SGETask instance created on the master node
- run SGETask.work()

On completion, SGETask on the master node will detect that
the job has left the queue, delete the temporary folder, and
return from SGETask.run()
"""

import os
import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle
import logging
import tarfile
import time


def _do_work_on_compute_node(work_dir, tarball=True, arrayjob=False):

    if tarball:
        # Extract the necessary dependencies
        # This can create a lot of I/O overhead when running many SGEJobTasks,
        # so is optional if the luigi project is accessible from the cluster node
        if(not arrayjob):
            _extract_packages_archive(work_dir)
        else:
            SGE_TASK_ID = int(os.environ['SGE_TASK_ID'])
            if(SGE_TASK_ID==1):
                _extract_packages_archive(work_dir)
                open(work_dir+'/packages_extracted.touch', 'a').close()
            else:
                timeout=60 #s to wait for SGE_TASK_ID 1 to untar deps before crashing
                while(not os.path.isfile(work_dir+'/packages_extracted.touch')
                      and timeout>0):
                    time.sleep(1)
                    timeout-=1

                if os.path.isfile(work_dir+'/packages_extracted.touch'):
                    if work_dir not in sys.path:
                        sys.path.insert(0, work_dir) #make sure we look in the untarred folders
                else:
                    JOB_ID = int(os.environ['JOB_ID'])
                    raise(Exception("{JOB}.{TASK} "
                                    "timed out while waiting for {JOB}.1 "
                                    "to untar dependencies.".format(
                                        JOB=JOB_ID, TASK=SGE_TASK_ID)))
                    exit(1)


    # Open up the pickle file with the work to be done
    os.chdir(work_dir)
    with open("job-instance.pickle", "rb") as f:
        job = pickle.load(f)
        
    #record hostname
    print("SGE job: {}:{}".format(int(os.environ['JOB_ID']), int(os.environ['SGE_TASK_ID'])))
    print("hostname: ", os.environ['HOSTNAME'])

    # Do the work contained
    job.work()


def _extract_packages_archive(work_dir):
    package_file = os.path.join(work_dir, "packages.tar")
    if not os.path.exists(package_file):
        return

    curdir = os.path.abspath(os.curdir)

    os.chdir(work_dir)
    tar = tarfile.open(package_file)
    for tarinfo in tar:
        tar.extract(tarinfo)
    tar.close()
    if '' not in sys.path:
        sys.path.insert(0, '')

    os.chdir(curdir)


def main(args=sys.argv):
    """Run the work() method from the class instance in the file "job-instance.pickle".
    """
    try:
        tarball = "--no-tarball" not in args
        arrayjob = "--arrayjob" in args
        # Set up logging.
        logging.basicConfig(level=logging.WARN)
        work_dir = args[1]
        assert os.path.exists(work_dir), "First argument to sge_runner.py must be a directory that exists"
        project_dir = args[2]
        sys.path.append(project_dir)
        _do_work_on_compute_node(work_dir, tarball, arrayjob)
    except Exception as e:
        # Dump encoded data that we will try to fetch using mechanize
        print(e)
        raise


if __name__ == '__main__':
    main()
