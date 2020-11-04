#! /bin/bash
#source /path/to/patched/gromacs/bin/GMXRC
export GMXLIB=../../../data/mutff

python SGE_tasks/absFE/Workflow_aligned_Complete.py --mdrun mdrun --mdrun_double mdrun --mdrun_opts=" -ntmpi 1 -notunepme " --toppath . --mdppath ./mdp --workers 8
