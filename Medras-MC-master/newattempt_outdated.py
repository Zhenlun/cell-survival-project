

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

from repairanalysis import medrasrepair
from repairanalysis import sddparser
from itertools import groupby

'''
individual break set field:
    0 is break set (index that starts at 0)
    1 is number of breaks
    2 is residual
    3 is number of misrepairs (include large misrepairs) (not the original pair)
    4 is number of large misrepairs (larger than 3 MBP)
    5 is inter-chromosome misrepairs (subset of index 3)
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
add_index_post = [2, 3]
add_index_pre = [1]

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
    
    #print(repair_data)
    return repair_data

def record_viability(num_simulations, sddpath):
    repair_data = run_simulation(sddpath)
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
   
    if plot_histogram: # Plot a histogram of the viability records
        # Define the bins such that 0 is always on the left
        bins = [0, 1, 2]
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
    plot_histogram = False
    path = os.getcwd()
    sddpath = path
    folder_rate = []
    folder_variance = []
    folder_deviation = []
    repeats = 100
    '''allfoldername = ['_inner_electron', '_inter_electron', '_outer_electron', 
                             '_inner_proton', '_inter_proton', '_outer_proton', 
                             '_inner_alpha', '_inter_alpha', '_outer_alpha']
    allfoldername = ['gy']'''
    #allfoldername =['proton']
    allfoldername = [str(i) for i in range(1, 101)]
    excludefoldername = ["proton"]
    #excludefoldername = ["gasdfg"]
    
    with open("cluster_single.txt", "w") as output_file:
        for root, subfolders, files in os.walk(sddpath):
            for foldername in subfolders:
                if any(sub in foldername for sub in allfoldername) and any(sub not in foldername for sub in excludefoldername):
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
                                    header, events = sddparser.parseSDDFile(file_path)
                                    folder_dose.append(header['Dose or Fluence'][1])
                    if len(folder_dose) > 0 and all_equal(folder_dose):
                        dose = folder_dose[0]
                    elif len(folder_dose) > 0 and not all_equal(folder_dose):
                        dose = np.mean(folder_dose)
                        #print('note that the folder contains files of different doses', file=output_file, flush=True)
                    else:
                        print('there are no files in this folder!', file=output_file, flush=True)
                    print([rate, std_dev, dose, folder_path], file=output_file, flush=True)
    
    
    '''for root, subfolders, files in os.walk(sddpath):
        for foldername in subfolders:
            if any(sub in foldername for sub in allfoldername) and any(sub not in foldername for sub in excludefoldername):
                folder_path = os.path.join(root, foldername)
                rate, var, std_dev = calculate_statistics(repeats, folder_path)
                folder_rate.append(rate)
                folder_variance.append(var)
                folder_deviation.append(std_dev)
                
                folder_dose = []
                for froot, fsubfolders, ffiles in os.walk(folder_path):
                        for file in ffiles:
                            #print(os.path.join(folder_path, file))
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
                print([rate, std_dev, dose, folder_path])'''



