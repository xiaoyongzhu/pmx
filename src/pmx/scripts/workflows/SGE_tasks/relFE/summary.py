#!/usr/bin/env python

import copy
import luigi
import os
import numpy as np
#import matplotlib as plt
from luigi.parameter import ParameterVisibility
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask #tuned for the owl cluster
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.analysis import Task_PL_analysis_aligned,Task_PL_analysis_aligned2crystal
from pmx.scripts.workflows.SGE_tasks.absFE.LinW.analysis import Task_WL_analysis


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_summary_aligned(SGETunedJobTask):
    #run on the login node
    run_locally = luigi.BoolParameter(
        default = True,
        significant=False,
        parsing=luigi.BoolParameter.EXPLICIT_PARSING,
        description="run locally instead of on the cluster")

    #Parameters:
    hosts = luigi.ListParameter(description='list of protein names to evaluate')
    ligands = luigi.ListParameter(description='list of ligand names to evaluate')
    n_repeats = luigi.IntParameter(description='Number of repeats',
                                   default=3)
    n_sampling_sims = luigi.IntParameter(description='Number of sampling simulations',
                                         default=1)

    #TODO: add default
    study_settings = luigi.DictParameter(significant=False,
        visibility=ParameterVisibility.HIDDEN,
        description='Dict of study stettings '
        'used to propagate settings to dependencies')

    #request 1 core
    n_cpu = luigi.IntParameter(default=1, significant=False,
                               visibility=ParameterVisibility.HIDDEN)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.WL_settings=copy.deepcopy(self.study_settings.get_wrapped())
        self.WL_settings['TIstates']=self.WL_settings['states']
        self.WL_settings['n_repeats']=self.n_repeats
        self.WL_settings['n_sampling_sims']=self.n_sampling_sims
        self.PL_settings=copy.deepcopy(self.study_settings.get_wrapped())
        self.PL_settings['n_repeats']=self.n_repeats
        self.PL_settings['n_sampling_sims']=self.n_sampling_sims
        #self.PL_settings=self.study_settings

        self.base_path = self.study_settings['base_path']
        self.outname="summary_aligned.txt"
        self.restrname="out_dg_{i}.dat"

        self.anafolderfmt_P="/analysis/repeat{i}"
        self.anafolderfmt_W="/analysis/repeat{i}"

    def work(self):

        def read_results(inP=False):
            rs=np.ndarray(self.n_repeats)
            for i in range(self.n_repeats):
                if(inP):
                    ana_folder=folder_path+self.anafolderfmt_P.format(i=i)
                else:
                    ana_folder=folder_path+self.anafolderfmt_W.format(i=i)
                with open(ana_folder+"/results.txt", 'r') as f:
                    for line in f:
                        if "BAR: dG" in line:
                            s = line.split()
                            rs[i]=float(s[3])
                            break

                if(inP): #analytical correction for this repeat
                    with open(folder_path+"/"+self.restrname.format(i=i), 'r') as f:
                        for line in f:
                            if("Restraint contribution to free energy" in line and
                               "kJ/mol" in line):
                                s = line.split()
                                rs[i]+=float(s[-2])
                                break

            dGpart = np.mean(rs)
            std = np.std(rs)
            return([dGpart,std])

        #dG in water
        inws={}
        p="water"
        for l in self.ligands:
            key=l
            folder_path = self.base_path+'/'+p+'/lig_'+l
            inws.update({key:read_results(inP=False)})

        #dG in protein
        inps={}
        for p in self.hosts:
            for l in self.ligands:
                key=p+' '+l
                folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                inps.update({key:read_results(inP=True)})

        #read analytical corrections for restraints
        anacorrs={}
        for p in self.hosts:
            for l in self.ligands:
                key=p+' '+l
                folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                cors=np.zeros(self.n_repeats)
                for i in range(self.n_repeats):
                    with open(folder_path+"/"+self.restrname.format(i=i), 'r') as f:
                        for line in f:
                            if("Restraint contribution to free energy " in line and
                               "kJ/mol" in line):
                                s=line.split()
                                cors[i]=float(s[-2])
                anacorrs.update({key:[np.mean(cors),np.std(cors)]})


        #print summary table
        with open(self.outname, 'w') as sf:

            print("{:^20s} \t{:^20s}   {:^20s}   {:^20s}   {:^20s}".format(
                            "host guest","ddG (kJ/mol)","dG in prot+restr" ,
                            "dG in water", "ana. corr."))
            sf.write("{:^20s} \t{:^20s}   {:^20s}   {:^20s}   {:^20s}\n".format(
                            "host guest","ddG (kJ/mol)","dG in prot+restr" ,
                            "dG in water", "ana. corr."))
            for p in self.hosts:
                for l in self.ligands:
                    key=p+' '+l
                    ddG = inws[l][0] - inps[key][0] #water - (protein + restr corr.)
                    sigma = np.sqrt(inps[key][1]**2 + inws[l][1]**2) #standard dev.
                    print("{:<20s}:\t{:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}".format(
                        key, ddG, sigma,
                        inps[key][0], inps[key][1],
                        inws[l][0], inws[l][1],
                        anacorrs[key][0], anacorrs[key][1]) )
                    sf.write("{:<20s}:\t{:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}   {:>8.2f} +- {:<8.2f}\n".format(
                        key, ddG, sigma,
                        inps[key][0], inps[key][1],
                        inws[l][0], inws[l][1],
                        anacorrs[key][0], anacorrs[key][1]) )

    def output(self):
        files=[self.outname]
        return([luigi.LocalTarget(os.path.join(self.base_path, f)) for f in files])

    def requires(self):
        tasks=[]

        #Ligand in Water
        p="water"
        for l in self.ligands:
            folder_path = self.base_path+'/'+p+'/lig_'+l
            for sTI in self.WL_settings['states']: #uses equil states for TI
                for i in range(self.WL_settings['n_repeats']):
                    tasks.append(Task_WL_analysis(
                        l = l, i = i,
                        study_settings = self.WL_settings,
                        folder_path = folder_path,
                        parallel_env=self.parallel_env))

        #Ligand in Protein
        for p in self.hosts:
            for l in self.ligands:
                folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                for sTI in self.PL_settings['TIstates']:
                    for i in range(self.PL_settings['n_repeats']):
                        tasks.append(Task_PL_analysis_aligned(
                            p = p, l = l, i = i,
                            study_settings = self.PL_settings,
                            folder_path = folder_path,
                            parallel_env=self.parallel_env))

        return(tasks)

class Task_summary_aligned2crystal(Task_summary_aligned):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #overwrite values
        self.outname="summary_aligned2crystal.txt"
        self.restrname="out_dg_aligned2crystal_{i}.dat"
        self.anafolderfmt_P="/analysis/aligned2crystal_repeat{i}"

    def requires(self):
        tasks=[]

        #Ligand in Water
        p="water"
        for l in self.ligands:
            folder_path = self.base_path+'/'+p+'/lig_'+l
            for sTI in self.WL_settings['states']: #uses equil states for TI
                for i in range(self.WL_settings['n_repeats']):
                    tasks.append(Task_WL_analysis(
                        l = l, i = i,
                        study_settings = self.WL_settings,
                        folder_path = folder_path,
                        parallel_env=self.parallel_env))

        #Ligand in Protein
        for p in self.hosts:
            for l in self.ligands:
                folder_path = self.base_path+'/prot_'+p+'/lig_'+l
                for sTI in self.PL_settings['TIstates']:
                    for i in range(self.PL_settings['n_repeats']):
                        tasks.append(Task_PL_analysis_aligned2crystal(
                            p = p, l = l, i = i,
                            study_settings = self.PL_settings,
                            folder_path = folder_path,
                            parallel_env=self.parallel_env))

        return(tasks)