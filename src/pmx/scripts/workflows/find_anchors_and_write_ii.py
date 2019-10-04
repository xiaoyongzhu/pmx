import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import scipy as sp
import scipy.optimize
import scipy.stats
from scipy.special import erf

#import tqdm


import MDAnalysis as md
from MDAnalysis.lib.distances import *
from MDAnalysis.analysis.distances import *

################################################################################
def _calc_restraint_dg_w_gromacs_limits(FCs,eqs):
        """Calculates effect of restraints in gas phase using limits of gromacs harmonic restraints.

        Returns restraint FE correction in kJ/mol
        This is the case for harmonic restraint potentials in gromacs.
        dist       [0, inf]
        ang        [0, Pi]
        dih-dih_0  [-Pi, Pi]
        """

        V0 = 1.66            # standard volume in nm^3
        kT = 2.4777097058    # @ 298K

        #-1^6 = 1;
        #Boresch, ... and Karplus put their Is in here
        # and they work out to I=-2 given infinite limits of constrant variables.
        # So in their case there is a np.sqrt(2*Pi*RT)^6 instead of np.sqrt(0.5*Pi*RT)^6
        integral_multiplier = np.power(np.pi*kT/2.0, 3.0)

        I_dist = -1 - erf( np.sqrt(FCs[0]/(2.0*kT)) * eqs[0] )
        I_ang = [0.0, 0.0]
        I_dih = [0.0, 0.0, 0.0]
        I_ang[0] = erf( np.sqrt(FCs[1]/(2.0*kT)) * ((eqs[1]*np.pi/180.0)-np.pi) ) - \
                   erf( np.sqrt(FCs[1]/(2.0*kT)) * (eqs[1]*np.pi/180.0) )
        I_ang[1] = erf( np.sqrt(FCs[2]/(2.0*kT)) * ((eqs[2]*np.pi/180.0)-np.pi) ) - \
                   erf( np.sqrt(FCs[2]/(2.0*kT)) * (eqs[2]*np.pi/180.0) )
        for i in range(3):
            I_dih[i] = erf( np.sqrt(FCs[3+i]/(2.0*kT)) * (-np.pi) ) - \
                       erf( np.sqrt(FCs[3+i]/(2.0*kT)) * (np.pi) )

        forceConstants = np.sqrt(FCs[0]*FCs[1]*FCs[2]*FCs[3]*FCs[4]*FCs[5])
        partition_func = np.power(eqs[0],2.0)*np.sin(eqs[1]*np.pi/180.0)*np.sin(eqs[2]*np.pi/180.0) * \
                         integral_multiplier * I_dist * I_ang[0] * I_ang[1] * \
                         I_dih[0] * I_dih[1] * I_dih[2] / forceConstants

        dg = -kT * np.log(8.0*np.power(np.pi,2.0)*V0 / partition_func)
        return(dg)

################################################################################
def harmonic(x, k, m, c):
    return(0.5*k*(x-m)**2 + c);
################################################################################
def harmonic_dih(x, k, m, c):
    new_x=np.where(x-m>180,  x-2*180, x)
    new_x=np.where(x-m<-180, x+2*180, new_x)
    return(0.5*k*(new_x-m)**2 + c);
################################################################################
def wrap_ang(x,m):
    new_x=np.where(x-m>90,  x-2*90, x)
    new_x=np.where(x-m<-90, x+2*90, new_x)
    return(new_x);
################################################################################
def wrap_dih(x,m):
    new_x=np.where(x-m>180,  x-2*180, x)
    new_x=np.where(x-m<-180, x+2*180, new_x)
    return(new_x);
################################################################################
def find_distribs(localindeces, relcoords, bscale=1.0, plot=False):
    dist=[]
    angA=[]
    angB=[]
    dihA=[]
    dihB=[]
    dihC=[]
    for fr in range(relcoords.shape[0]):                        # go over each frame
        pos = relcoords[fr,localindeces,:]
        dist.append(calc_bonds(pos[0],pos[3]))
        angA.append(calc_angles(pos[1],pos[0],pos[3]))
        angB.append(calc_angles(pos[0],pos[3],pos[4]))
        dihA.append(calc_dihedrals(pos[2],pos[1],pos[0],pos[3]))
        dihB.append(calc_dihedrals(pos[1],pos[0],pos[3],pos[4]))
        dihC.append(calc_dihedrals(pos[0],pos[3],pos[4],pos[5]))

    data=[dist, angA, angB, dihA, dihB, dihC]
    for i in range(len(data)):
        data[i]=np.array(data[i])
        if(i==0):
            data[i]/=10.0 #nm to A
        else:
            data[i]*=180.0/np.pi #rad to deg

    dist_h = np.histogram(data[0], bins=int(200/bscale), range=(0,2.0), density=True)
    angA_h = np.histogram(data[1], bins=int(180/bscale), range=(0,180), density=True)
    angB_h = np.histogram(data[2], bins=int(180/bscale), range=(0,180), density=True)
    dihA_h = np.histogram(data[3], bins=int(360/bscale), range=(-180,180), density=True)
    dihB_h = np.histogram(data[4], bins=int(360/bscale), range=(-180,180), density=True)
    dihC_h = np.histogram(data[5], bins=int(360/bscale), range=(-180,180), density=True)



    #convert bins into x coordinate
    dists=[dist_h, angA_h, angB_h, dihA_h, dihB_h, dihC_h]
    for i in range(len(dists)):
        d=dists[i]
        dists[i] = [d[0], d[1][:-1]+(d[1][1]-d[1][0])*0.5]

        if(i>0):
            m=dists[i][1][np.argmax(dists[i][0])] #coord of max height of distribution
            if(i<3): #wrap angles
                data[i]=wrap_ang(data[i],m)
                dists[i][1]=wrap_ang(dists[i][1],m)
            else:    #wrap dihedrals
                data[i]=wrap_dih(data[i],m)
                dists[i][1]=wrap_dih(dists[i][1],m)
            #sort so distribution is in assending order of coordinate
            c = np.argsort(dists[i][1])
            dists[i][1]=dists[i][1][c]
            dists[i][0]=dists[i][0][c]

    #find mean and std (after wraping to PBC)
    equil_vals=[np.mean(d) for d in data]
    sigmas=[np.std(d) for d in data]


    dist_names=["distance", "angA", "angB", "dihA", "dihB", "dihC"]
    units=["nm", "deg", "deg", "deg", "deg", "deg"]
    units_FC=["kJ/(mol*nm^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)", "kJ/(mol*rad^2)"]

    kT = 2.4777097058
    pot = lambda x,k,m,c: 0.5*k*((x-m)*np.pi/180)**2+c
    PMF = lambda y: -kT*np.log(y)
    guess_FC=kT/(np.array(sigmas)**2) #in kJ/mol * ang^-2 and deg^-2

    if(plot):
        plt.figure(figsize=(6, 18))
        sp_grid = gridspec.GridSpec(6, 2)

    for i in range(6):
        pmf=PMF(dists[i][0])

        if(plot):
            plt.subplot(sp_grid[i,0])
            plt.title("%s distribution"%dist_names[i])
            plt.xlabel("%s (%s)"%(dist_names[i],units[i]))
            plt.ylabel("probability")
            plt.plot(dists[i][1], dists[i][0], '-b')
            #reference normal distribution based on mean and std
            rv=sp.stats.norm(loc=equil_vals[i], scale=sigmas[i]).pdf(dists[i][1])
            plt.plot(dists[i][1], rv, '--r')

            #PMF of data
            plt.subplot(sp_grid[i,1])
            plt.title("%s PMF"%dist_names[i])
            plt.xlabel("%s (%s)"%(dist_names[i],units[i]))
            plt.ylabel("U (kJ/mol)")
            plt.plot(dists[i][1], pmf, '-b', label="PMF")

        #sanitize infs and nans for fitting/plotting
        mask = np.isfinite(pmf)
        fit_y = pmf[mask]
        fit_x = dists[i][1][mask]

        #fit k
        popt=[guess_FC[i], equil_vals[i], -kT*np.log(1/(sigmas[i]*np.sqrt(2*np.pi)))]


        k=popt[0]
        if(i>0):
            k = k*((180.0/np.pi)**2);

        if plot:
            lbl="harmonic U\nk=%6.2f %s"%(k,units_FC[i])
            plt.plot(dists[i][1], harmonic(dists[i][1],*popt), '--r',\
                     label=lbl)
            plt.legend()

            #set plot y limits to still see the PMF
            spread=np.ptp(fit_y)
            plt.ylim((np.min(fit_y)-0.05*spread, np.max(fit_y)+0.4*spread))

        #prep k & m for output
        guess_FC[i]=k;
        equil_vals[i]=popt[1];

    if plot:
        plt.tight_layout()
        plt.savefig("restraint_coord_distrib.png", dpi=300)
        #plt.show()


    Dsum = 0
    for i in range(6):
        cdf=sp.stats.norm(loc=equil_vals[i], scale=sigmas[i]).cdf
        D,p=sp.stats.kstest(data[i],cdf,N=100)
        #print(D,p)
        Dsum+=D

    #if equilibrium angles near limits, penalize Dsum
    for i in range(len(equil_vals)):
        if(i>0):
            limits=[]
            if(i<3): #angle
                limits=[0.0,180.0]
            else:    #dihedral
                limits=[-180.0,180.0]
            if(np.any(np.greater(4*kT, harmonic(limits, guess_FC[i], equil_vals[i], 0.0)))):
                Dsum+=1.0 #if potential at limits is < 4 kT -> bad
                if(plot):
                    print("Penalizing Dsum for a low barier at limit!")

    if(plot):
        return(guess_FC,equil_vals);
    return(Dsum)

################################################################################
def find_restraints(structA="dumpA.gro", structB="dumpB.gro",
                    ref="averageA.gro",
                    trajA="all_eqA_fit.xtc", trajB="all_eqB_fit.xtc",
                    out="ii.itp", an_cor_file="out_dg.dat",
                    skip=1, log=False):
    
    
    #load trajectories
    A = md.Universe(structA,trajA)
    #lig_cand = A.select_atoms("resnum 2 and not (name H*)")
    lig_cand = A.select_atoms("resname MOL and not (name H*)")

    lcog=lig_cand.centroid()

    ref = md.Universe(ref)
    #lig_ref = ref.select_atoms("resnum 2 and not (name H*)")
    lig_ref = ref.select_atoms("resname MOL and not (name H*)")

    B = md.Universe(structB,trajB)
    #prot_cand = B.select_atoms("resnum 1 and not (name H*) and point %f %f %f %f"%(lcog[0], lcog[1], lcog[2], 15.0))
    prot_cand = B.select_atoms("protein and not (name H*) and point %f %f %f %f"%(lcog[0], lcog[1], lcog[2], 15.0))
    

    #ligand variance in A
    displ_mat_A=[]
    for frame in A.trajectory:                        # go over each frame in state A
        if(frame.frame%skip==0):
            lig_pos = lig_cand.positions
            lig_ref_pos = lig_ref.positions
            displ_mat_A.append(lig_pos-lig_ref_pos)
    displ_mat_A=np.array(displ_mat_A)
    var_A = np.var(displ_mat_A, axis=0)

    displ_mat_B=[]
    for frame in B.trajectory:                        # go over each frame in state B
        if(frame.frame%skip==0):
            prot_pos = prot_cand.positions
            lig_ref_pos = lig_ref.positions
            displ_mat_B.append(prot_pos[:,np.newaxis,:]-lig_ref_pos[np.newaxis,:,:])

    #ligand variance in B
    displ_mat_B=np.array(displ_mat_B)
    var_B = np.var(displ_mat_B, axis=0)
    var_combined = var_B + var_A[np.newaxis,:,:]

    #total variance
    varsq = np.linalg.norm(var_combined, axis=2)

    #sort by variance
    flat_varsq = varsq.reshape((-1))
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
    parser = argparse.ArgumentParser(description='Finds best anchor atoms, fits the force constants and writes the restraints itp file.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-sA', dest='structA', type=str, default="dumpA.gro",
                       help='input structure for A')
    parser.add_argument('-sB', dest='structB', type=str, default="dumpB.gro",
                       help='input structure for B')
    parser.add_argument('-r', dest='ref', type=str, default="averageA.gro",
                       help='input reference average structure of A')
    parser.add_argument('-fA', dest='trajA', type=str, default="all_eqA_fit.xtc",
                       help='input trajectory')
    parser.add_argument('-fB', dest='trajB', type=str, default="all_eqB_fit.xtc",
                       help='input trajectory')
    parser.add_argument('--oii', dest='out', type=str, default="ii.itp",
                       help='output restraint topology')
    parser.add_argument('-dg', dest='an_cor_file', type=str, default="out_dg.dat",
                       help='output restraint topology')
    parser.add_argument('--skip', dest='skip', type=int, default=1,
                       help='analyse every this many frames')

    args = parser.parse_args()

    find_restraints(args.structA, args.structB, args.ref,
                    args.trajA, args.trajB, args.out, args.an_cor_file,
                    args.skip, log=True)


    exit(0);

