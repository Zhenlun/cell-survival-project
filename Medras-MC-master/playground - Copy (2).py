
# ############################################################################
# 
# This software is made freely available in accordance with the simplifed BSD
# license:
# 
# Copyright (c) <2017>, <Stephen McMahon>
# All rights reserved
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contacts: Stephen McMahon,	stephen.mcmahon@qub.ac.uk
# 
# ############################################################################
#
# References
#
# 1. McMahon et al, Scientific Reports, 2016, 6, 33290
# 2. McMahon et al, Scientific Reports, 2017, 7, 10790
#
# ############################################################################


import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns

import os
import numpy as np
import scipy as sp
import math
from scipy.integrate import odeint
from scipy.stats import poisson
import sys



from repairanalysis import medrasrepair
from repairanalysis import sddparser
from itertools import groupby

path = os.getcwd()
#sddpath = path


#print('Fidelity analysis')
#fide = medrasrepair.repairSimulation(sddpath,'Fidelity', 24)
#fide = fide[0]
#fide = fide.strip().split()

#print('\n\nMisrepair spectrum analysis')
#medrasrepair.repairSimulation(sddpath,'Spectrum')



'''
individual break set field:
    0 is file name
    1 is break set
    2 is break count
    3 is misrepair (in DSB? in fractional?)
    4 is misrepair standard deviation
    5 is inter-chromosome rate
    
    print('File\tBreak Set\tBreak Count\tMisrepair\tStdev\tInter-Chromosome Rate', end='')
    

summary field:
    0 is file name
    1 is break set (summary)
    2 is total break count
    3 is complexity (unit? how is it counted?)
    4 is average break count
    5 is average break count's standard deviation
    6 is average misrepairs
    7 is misrepair standard deviation
    8 is file average (of what? error rate? whatever that is)
    
	smry = (fileName+'\tSummary\t'+str(totalBreaks)+'\t'+str(complexity)+'\t'+str(averageBreaks)+
		   '\t'+str(breakStdev)+'\t'+str(averageMisrep)+'\t'+str(misrepStdev)+'\t'+
		   str(fileAverages[0]) )
'''

#print('\n\nSeparation analysis')
#sepa = medrasrepair.repairSimulation(sddpath,'Separation')


#sys.stdout = open("medras_spectrum_analysis.output", 'w')
#medrasrepair.repairSimulation(sddpath,'Spectrum')
#sys.stdout.close()
#print('\n\nMisrepair spectrum analysis')
#spec = medrasrepair.repairSimulation(sddpath,'Spectrum', 24)
#spec = spec[0]
'''
individual break set field:
    0 is break set (index that starts at 0)
    1 is number of breaks
    2 is residual (??? guess is number of DSB that remains)
    3 is number of misrepairs (include large misrepairs?)
    4 is number of large misrepairs (larger than 3 MBP)
    5 is inter-chromosome misrepairs
    6 is the # of dicentrics
    7 is the # of rings
    8 is the # of excess linear fragments
    9 is the # of total aberrations
    10 is the viability

'''







defaultSurvParams= {'apoptoticRate': 0.01117, 'mitoticRate': 0.0141, 'baseRate': 0.000739} 
defaultDNAParams = {'sigma': 0.04187, 'NHEJFidelity': 0.98537,'MMEJFidelity': 0.4393, 
    						'fastRepair': 2.081,'fastFociDelay': 8.079, 'slowRepair': 0.2604, 
    						'slowFociDelay': 0.4053,'verySlowRepair': 0.008462, 
    						'complexFrac': 0.4337, 'pointMutationRate': 0.04605,   'failFrac': 0.7364,
    						'rNuc': 4.282 }
defaultCell = {'dna':6100.0,'chromosomes':46,'repair':0,'G1Arrest':1,'phase':0,'gene':0}  


add_index_aber = [9]
#add_index_aber = [10, 11, 14, 15]
add_index_post = [2, 3]
add_index_pre = [1]


plot_survival_rate = False
plot_phase_rate = False

functional_g1arrest = True
cell_phase = [1,2,3]
LET_value = 0
dose_value = 2
particle = 2112 #PDG code

currentstdout = sys.__stdout__

def run_simulation(sddp):
    #print('\n\nMisrepair spectrum analysis')
    
    sys.stdout = open("medras_spectrum_analysis.output", 'w')
    medrasrepair.repairSimulation(sddp,'Spectrum')
    sys.stdout.flush()
    sys.stdout = sys.__stdout__
    
    with open("medras_spectrum_analysis.output", 'r') as spectrum_file:
        repair_data = [line.split() for line in spectrum_file if line[0].isdigit()]
    
    return repair_data



def g1_survival(aber, prebreak, postbreak):
    """
    Calculate G1 survival rate based on aberrations and breaks.
    """
    apoptotic_survival = np.exp(-postbreak * defaultSurvParams['apoptoticRate']) if functional_g1arrest else 1.0
    aberration_survival = np.exp(-aber)
    base_survival = np.exp(-prebreak * defaultSurvParams['baseRate'])
    
    return apoptotic_survival * aberration_survival * base_survival
        
    
    
def g2_survival(aber, postbreak):
    """
    Calculate G2 survival rate based on aberrations and post-repair breaks.
    """
    mitotic_break = min(postbreak, 20)
    mitotic_survival = np.exp(-mitotic_break * defaultSurvParams['mitoticRate'])
    aberration_survival = np.exp(-aber)
    
    return mitotic_survival * aberration_survival
    
    
    
def m_survival(aber, prebreak, postbreak):
    """
    Calculate M survival rate based on aberrations and breaks.
    """
    mitotic_break = postbreak
    mitotic_survival = np.exp(-mitotic_break * defaultSurvParams['mitoticRate'])
    g1_survival_rate = g1_survival(aber, prebreak, postbreak)
    
    return mitotic_survival * g1_survival_rate



def process_spectrum(spectrum, phase):
    return_rate = []
    
    for spec in spectrum:
        # Assuming add_index_aber, add_index_post, and add_index_pre are defined elsewhere
        aberration = sum(int(spec[i]) for i in add_index_aber if len(spec) > 7)
        postbreak = sum(int(spec[i]) for i in add_index_post if len(spec) > 7)
        prebreak = sum(int(spec[i]) for i in add_index_pre if len(spec) > 7)
        
        if phase == 1:
            rate = g1_survival(aberration, prebreak, postbreak)
        elif phase == 2:
            rate = g2_survival(aberration, postbreak)
        elif phase == 3:
            rate = m_survival(aberration, prebreak, postbreak)
        else:
            raise ValueError("Please enter a valid phase!")
        
        return_rate.append(rate)
    
    return return_rate
    
        
    
def get_survival_rates(sddp):
    spectrum = run_simulation(sddp)
    
    all_phase_rate = [process_spectrum(spectrum, phase) for phase in cell_phase]
    return all_phase_rate




def go_through_folder(subfolder):
    sddpath = os.path.join(os.getcwd(), subfolder)
    result = get_survival_rates(sddpath)
    return result

def plot_results(results):
    def linq1(x):
        a = 2.01 
        b = 0.011
        return np.exp(-a*x-b*x*x)

    def linq2(x):
        a = 0.5
        b = 0.064
        return np.exp(-a*x-b*x*x)

    xax = np.linspace(0, 4)
    plt.yscale("log")
    plt.plot(xax, linq1(xax), label = "chaudhary14, AG01522 cells")
    plt.plot(xax, linq2(xax), label = "chaudhary14, U-87 cells")
    for i, result in enumerate(results):
        plt.errorbar(0.094, result[0][0], yerr=result[0][1], fmt='r.')#, label=f"this work {i+1}")
    plt.legend()
    plt.title("Model Survival Rate In Comparison To PIDE Data")
    plt.xlabel("Dose (Gy)")
    plt.ylabel("Survival Rate")
    plt.show()
    
    
def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)
    
if __name__ == "__main__":
    sddpath = os.getcwd()
    repeats = 10000
    
    for root, subfolders, files in os.walk(sddpath):
        for foldername in subfolders:
            if any(sub in foldername for sub in ['_inner_electron', '_inter_electron', '_outer_electron', 
                                     '_inner_proton', '_inter_proton', '_outer_proton', 
                                     '_inner_alpha', '_inter_alpha', '_outer_alpha']):

                g1_result = []
                g2_result = []
                m_result = []
                subfolder = os.path.join(root, foldername)
                splength = len(run_simulation(subfolder))
                if splength > 0:
                    sim_to_run = math.ceil(repeats / splength)
                    for i in range(sim_to_run):
                        results = go_through_folder(subfolder)
                        for i in range(len(results[0])):
                            g1_result.append(results[0][i])
                            g2_result.append(results[1][i])
                            m_result.append(results[2][i])
                else:
                    print(splength)
                g1_rate = np.mean(g1_result)
                g1_std = np.std(g1_result)
                g2_rate = np.mean(g2_result)
                g2_std = np.std(g2_result)
                m_rate = np.mean(m_result)
                m_std = np.std(m_result)
                
                folder_path = subfolder
                folder_dose = []
                for froot, fsubfolders, ffiles in os.walk(folder_path):
                    for file in ffiles:
                        if file.endswith('.txt'):
                            file_path = os.path.join(folder_path, file)
                            #print(file_path)
                            header, events = sddparser.parseSDDFile(file_path)
                            folder_dose.append(header['Dose or Fluence'][1])
                if len(folder_dose) > 0 and all_equal(folder_dose):
                    dose = folder_dose[0]
                elif len(folder_dose) > 0 and not all_equal(folder_dose):
                    dose = np.mean(folder_dose)
                    
                print([g1_rate, g1_std, dose, subfolder], flush = True)


'''
result = statistics(repeats, path)
#result = [[0.10279999999999999, 0.027461973709112752], [0.1013, 0.02540295258429618], [0.1029, 0.02542813402513051]]
sys.stdout.flush()
sys.stdout = sys.__stdout__
print(result)


def linq1(x):
    a = 2.01 
    b = 0.011
    return np.exp(-a*x-b*x*x)

def linq2(x):
    a = 0.5
    b = 0.064
    return np.exp(-a*x-b*x*x)


plt.figure()
xax = np.linspace(0, 4)
plt.yscale("log")
plt.plot(xax, linq1(xax), label = "chaudhary14, AG01522 cells")
plt.plot(xax, linq2(xax), label = "chaudhary14, U-87 cells")
plt.errorbar(0.094, result[0][0],yerr = result[0][1], fmt='r.', label = "this work")
plt.legend()
plt.title("Model Survival Rate In Comparison To PIDE Data")
plt.xlabel("Dose (Gy)")
plt.ylabel("Survival Rate")
plt.show()'''
