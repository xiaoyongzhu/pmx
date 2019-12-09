import argparse
import numpy as np
import MDAnalysis as md
from MDAnalysis.lib.distances import *
from MDAnalysis.analysis.distances import *
from pmx.scripts.workflows.find_anchors_and_write_ii import find_distribs, _calc_restraint_dg_w_gromacs_limits

################################################################################
def find_restraints_align2crystal(struct="dumpD.gro",
                    traj="all_eqD_fit.xtc",
                    out="ii.itp", an_cor_file="out_dg.dat",
                    skip=1, log=False):


    #load trajectories
    A = md.Universe(struct,traj)
    #lig_cand = A.select_atoms("resnum 2 and not (name H*)")
    lig_cand = A.select_atoms("resname MOL and not (name H*)")

    lcog=lig_cand.centroid()

    prot_cand = A.select_atoms("backbone and not (name H*) and point %f %f %f %f"%(lcog[0], lcog[1], lcog[2], 15.0))


    #ligand variance in A
    displ_mat=[]
    for frame in A.trajectory:                        # go over each frame in state A
        if(frame.frame%skip==0):
            lig_pos = lig_cand.positions
            prot_pos = prot_cand.positions
            displ_mat.append(prot_pos[:,np.newaxis,:]-lig_pos[np.newaxis,:,:])

    #ligand variance in B
    displ_mat=np.array(displ_mat)
    var = np.var(displ_mat, axis=0)

    #total variance
    varsq = np.linalg.norm(var, axis=2)

    #sort by variance
    #flat_varsq = varsq.reshape((-1))
    s = np.argsort(varsq,  axis=None)

    #filter best candidates
    L=[]
    P=[]
    i=0
    lim=6
    while (len(L)<lim or len(P)<lim) and i<len(s):
        p=int(s[i]/varsq.shape[1])
        l=int(s[i]-p*varsq.shape[1])
        #ppos=prot_cand[p].position
        #lpos=lig_ref[l].position
        #if(np.linalg.norm(ppos-lpos)>5.0**2): #A
        #    continue; #pair too far, will produce very strong angle/dih force constants

        p=prot_cand[p].index
        l=lig_cand[l].index
        if(p not in P and len(P)<lim):
            P.append(p)
        if(l not in L and len(L)<lim):
            L.append(l)
        i+=1

    #load the trajectory of relevant atom positions in A to memory so we don't have to read trajectory many times
    relevant = md.core.groups.AtomGroup(P+L, A)
    relevantpos=np.zeros((int(len(A.trajectory)/skip), len(relevant), 3))
    i=0
    for frame in A.trajectory:                        # go over each frame
        if(frame.frame%skip==0):
            relevantpos[i,:,:]=relevant.positions
            i+=1

    #loop through all permutations of valid candidates
    possible_anchors=[]
    for p1 in range(len(P)):
        for p2 in range(p1,len(P)):
            if(p2==p1):
                continue
            for p3 in range(p2,len(P)):
                if(p3==p1 or p3==p2):
                    continue
                for l1 in range(len(L)):
                    for l2 in range(l1,len(L)):
                        if(l2==l1):
                            continue
                        for l3 in range(l2,len(L)):
                            if(l3==l1 or l3==l2):
                                continue
                            possible_anchors.append([p1,p2,p3,len(P)+l1,len(P)+l2,len(P)+l3])

    if(not possible_anchors):
        raise ValueError("No potential anchors found.")

    #find the most gaussian distributions among the permutations
    last_anchors=[]
    last_Dsum=1e10
    #for i in tqdm.tqdm(range(len(possible_anchors))):
    for i in range(len(possible_anchors)):
        test_anchors=possible_anchors[i]
        Dsum=find_distribs(test_anchors, relevantpos, plot=False)
        if(Dsum<last_Dsum):
            last_anchors=test_anchors
            last_Dsum=Dsum

    #output results and plot fits
    final_anchors=[relevant[i].index for i in last_anchors]
    if(log):
        print("Final Dsum:%.4f\t final anchors"%last_Dsum, final_anchors)
    FCs,eqs=find_distribs(last_anchors, relevantpos, plot=True)

    #add 1 to all anchor indices, as gmx indexing starts at 1, not 0.
    final_anchors=[relevant[i].index + 1 for i in last_anchors]

    #write restraint topology
    with open(out, "w") as ii:
        print("", file=ii)
        print(" [ intermolecular_interactions ]", file=ii)
        print(" [ bonds ]", file=ii)
        print("%6d %6d %6d %14.6f %14.6f %14.6f %14.6f"%(\
                final_anchors[0],final_anchors[3],\
                6, eqs[0],0.0,eqs[0],FCs[0]), file=ii)
        print("", file=ii)
        print(" [ angles ]", file=ii)
        print("%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f"%(\
                final_anchors[1],final_anchors[0],final_anchors[3],\
                1, eqs[1],0.0,eqs[1],FCs[1]), file=ii)
        print("%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f"%(\
                final_anchors[0],final_anchors[3],final_anchors[4],\
                1, eqs[2],0.0,eqs[2],FCs[2]), file=ii)
        print("", file=ii)
        print(" [ dihedrals ]", file=ii)
        print("%6d %6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f"%(\
                final_anchors[2],final_anchors[1],final_anchors[0],final_anchors[3],\
                2, eqs[3],0.0,eqs[3],FCs[3]), file=ii)
        print("%6d %6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f"%(\
                final_anchors[1],final_anchors[0],final_anchors[3],final_anchors[4],\
                2, eqs[4],0.0,eqs[4],FCs[4]), file=ii)
        print("%6d %6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f"%(\
                final_anchors[0],final_anchors[3],final_anchors[4],final_anchors[5],\
                2, eqs[5],0.0,eqs[5],FCs[5]), file=ii)

    #write analytical corection to FE because of restraints
    dg=_calc_restraint_dg_w_gromacs_limits(FCs,eqs)
    with open(an_cor_file, "w") as dgf:
        print('Restraint contribution to free energy (w gmx limits): %3.4f kJ/mol' % dg, file=dgf)
        print('Restraint contribution to free energy (w gmx limits): %3.4f kcal/mol' % (dg/4.184), file=dgf)
    if(log):
        print('Restraint contribution to free energy (w gmx limits): %3.4f kJ/mol' % dg)
        print('Restraint contribution to free energy (w gmx limits): %3.4f kcal/mol' % (dg/4.184))


################################################################################
#start of execution
if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Finds best anchor atoms, '
                                     'fits the force constants and writes the '
                                     'restraints itp file based on an aligned '
                                     'trajectory.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', dest='struct', type=str, default="dumpA.gro",
                       help='input structure for A')
    parser.add_argument('-f', dest='traj', type=str, default="all_eqA_fit.xtc",
                       help='input trajectory')
    parser.add_argument('--oii', dest='out', type=str, default="ii.itp",
                       help='output restraint topology')
    parser.add_argument('-dg', dest='an_cor_file', type=str, default="out_dg.dat",
                       help='output restraint topology')
    parser.add_argument('--skip', dest='skip', type=int, default=1,
                       help='analyse every this many frames')

    args = parser.parse_args()

    find_restraints_align2crystal(args.struct, args.traj,
                                  args.out, args.an_cor_file,
                                  args.skip, log=True)


    exit(0);

