
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
import math
import numpy as np
import scipy as sp
from scipy.integrate import odeint
from scipy.stats import poisson
import sys
import io
from contextlib import redirect_stdout

from concurrent.futures import ThreadPoolExecutor

from repairanalysis import medrasrepair
from repairanalysis import sddparser
from itertools import groupby


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

'''def run_simulation(sddpath):
    # Misrepair spectrum analysis
    sys.stdout = open("medras_spectrum_analysis.output", 'w')
    medrasrepair.repairSimulation(sddpath, 'Spectrum')
    sys.stdout.flush()
    sys.stdout = sys.__stdout__
    
    with open("medras_spectrum_analysis.output", 'r') as spectrum_file:
        repair_data = [line.split() for line in spectrum_file if line[0].isdigit()]
    
    return repair_data'''


def run_simulation(sddp):
    # Create an in-memory text stream
    memory_file = io.StringIO()
    
    # Redirect stdout to the in-memory stream
    with redirect_stdout(memory_file):
        medrasrepair.repairSimulation(sddp, 'Spectrum')
    
    # Reset the stream position to the beginning
    memory_file.seek(0)
    
    # Read the data from the in-memory stream
    repair_data = [line.split() for line in memory_file if line[0].isdigit()]
    
    # Close the in-memory stream
    memory_file.close()
    
    return repair_data



def record_viability(num_simulations, sddpath):
    repair_data = run_simulation(sddpath)
    #print(repair_data)
    #print(len(repair_data))
    viability_records = []
    if len(repair_data) > 0:
        # Run the simulation num_simulations/len(repair_data) times, rounded up
        simulations_to_run = math.ceil(num_simulations / len(repair_data))
        for _ in range(simulations_to_run):
            repair_data = run_simulation(sddpath)
            for data_set in repair_data:
                # Check if the last element is "misrepair!", indicating a viability of 1
                if data_set[-1] == "misrepair!":
                    viability_records.append('1')
                else:
                    # Otherwise, record the actual last element as the viability
                    viability_records.append(data_set[-1])
    else:
        return 'check folder!!!'
    return np.array(viability_records, dtype=int)

def plot_histogram_and_calculate_survival_rate(viability_records):
    # Define the bins such that 0 is always on the left
    bins = [0, 1, 2]
    
    # Plot a histogram of the viability records
    plt.figure()
    plt.hist(viability_records, bins=bins, align='left', rwidth=0.8, color='blue', alpha=0.7)
    plt.title('Histogram of Viability Outcomes')
    plt.xticks([0, 1])
    plt.xlabel('Viability')
    plt.ylabel('Frequency')
    plt.show()
    
    # Calculate the survival rate
    survival_rate = np.count_nonzero(viability_records == 1) / len(viability_records)
    return survival_rate

def calculate_survival_rate(num_simulations, sddpath):
    viability_results = record_viability(num_simulations, sddpath)
    #print(len(viability_results))
    survival_rate = plot_histogram_and_calculate_survival_rate(viability_results)
    return survival_rate

def calculate_statistics(number_of_trials, sddpath):
    survival_rate = calculate_survival_rate(number_of_trials, sddpath)
    n = number_of_trials
    p = survival_rate
    variance = p * (1 - p) / n
    standard_deviation = np.sqrt(variance)
    return survival_rate, variance, standard_deviation

def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


        
if __name__ == "__main__":
    path = os.getcwd()
    sddpath = path
    folder_rate = []
    folder_variance = []
    folder_deviation = []
    repeats = 10000

    for root, subfolders, files in os.walk(sddpath):
        for foldername in subfolders:
            if any(sub in foldername for sub in ['_inner_electron', '_inter_electron', '_outer_electron', 
                                     '_inner_proton', '_inter_proton', '_outer_proton', 
                                     '_inner_alpha', '_inter_alpha', '_outer_alpha']):
                folder_path = os.path.join(root, foldername)
                rate, var, std_dev = calculate_statistics(repeats, folder_path)
                folder_rate.append(rate)
                folder_variance.append(var)
                folder_deviation.append(std_dev)
                
                
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
                    print('note that the folder contains files of different doses')
                else:
                    print('there are no files in this folder!')
                print([rate, std_dev, dose, folder_path])



