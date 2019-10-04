import argparse
import numpy as np
import MDAnalysis as md
from MDAnalysis.lib.distances import *
from MDAnalysis.analysis.distances import *



def find_avg_struct(struct, traj, out, log=False):
    A = md.Universe(struct, traj)                 # always start with a Universe
    a = A.select_atoms("all")

    avg_mat=np.zeros_like(a.positions)
    nframes=0
    for frame in A.trajectory:                        # go over each frame
        avg_mat += a.positions
        nframes+=1
    avg_mat/=nframes

    ref = md.Universe(struct)
    ref.atoms.positions=avg_mat

    with md.Writer(out, multiframe=False, bonds=None, n_atoms=a.n_atoms) as OUT:
        OUT.write(ref)
        if(log):
            print("Wrote average structure into "+out)
    



if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Finds the average structure of a trajectory.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', dest='struct', type=str, required=True,
                       help='input structure')
    parser.add_argument('-f', dest='traj', type=str, required=True,
                       help='input trajectory')
    parser.add_argument('-o', dest='out', type=str, default="avg.gro",
                       help='output structure')

    args = parser.parse_args()

    find_avg_struct(args.struct, args.traj, args.out, log=True)

    exit(0);
