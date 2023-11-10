
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
from scipy.integrate import odeint
from scipy.stats import poisson
import sys



from repairanalysis import medrasrepair
from repairanalysis import plotAberrations


path = os.getcwd()
sddpath = path


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


repeats = 100
plot_survival_rate = False
plot_phase_rate = False

functional_g1arrest = True
cell_phase = [1,2,3]
LET_value = 0
dose_value = 2
particle = 2112 #PDG code



def run_simulation():
    #print('\n\nMisrepair spectrum analysis')
    
    sys.stdout = open("medras_spectrum_analysis.output", 'w')
    medrasrepair.repairSimulation(sddpath,'Spectrum')
    sys.stdout.flush()
    
    spectrum_file = open("medras_spectrum_analysis.output", 'r')
    lines = spectrum_file.readlines()
    repair_data = []
    
    for line in lines:
        if line[0].isdigit():
            repair_data.append(line.split())
            
    spectrum_file.close()
    #os.remove("medras_spectrum_analysis.output")
    
    #print(repair_data)
    return repair_data


def calculatesurvival(phase):
 
    spectrum = run_simulation()
    return_rate = []
    
    test_rate = 0
    for index in range(len(spectrum)):
        
        #for each individual set, will want to run the simulation multiple times
        #and obtain an average and error for each value
        
        
        #fidelity = fide[index]
        spec = spectrum[index]
        
        
        #lethal aberration
        aberration = 0
        for i in add_index_aber:
            if len(spec)>7:
                aberration += int(spec[i])

        
        #post repair breaks
        postbreak = 0
        for i in add_index_post:
            if len(spec)>7:
                postbreak += int(spec[i])
    
        
        #pre repair breaks
        prebreak = 0
        for i in add_index_pre:
            if len(spec)>7:    
                prebreak += int(spec[i])
            
        
        if int(phase) == 1:
            rate = g1_survival(aberration, prebreak, postbreak)
            #print("g1 survival for index " + str(index) +": " + str(rate))
            return_rate.append(rate)
            
        elif phase == 2:
            rate = g2_survival(aberration, postbreak)
            #print("g2 survival for index " + str(index) +": " + str(rate))
            return_rate.append(rate)
        
        elif phase == 3:
            rate = m_survival(aberration, prebreak, postbreak)
            #print("mitosis survival for index " + str(index) +": " + str(rate))
            return_rate.append(rate)
            
        else:
            print("please enter a valid phase!")
            
    if return_rate == []:
            
        sys.stdout.write(str(phase == 1))
        sys.stdout.flush()
            
        sys.stdout.write(str(test_rate))
        sys.stdout.flush()
            
        raise ValueError
            
    #print(return_rate)
    return return_rate
    
        
   
#returnConditions.append([currBreaks, misrepairedDSB, lethalAberrations, visibleAberrations, totalMut, self.totalDSBs])


def g1_survival(aber, prebreak, postbreak):
    '''
    how to figure out if g1 arrest is functional? currently no way to tell if the cell has functional g1arrest
						    
                            break at proliferation = postrepair breaks
                            apoptoticRate = full apoptotic rate
                            erate = base apoptotic rate
                            code for lethal aberration
						    totalbreaks = prerepair breaks 
							baserate
    '''
	
	
	
	
    if functional_g1arrest:
        apoptoticsurvival = np.exp(-(postbreak * defaultSurvParams['apoptoticRate']))
        #print(apoptoticsurvival)
    else:
        apoptoticsurvival = 1.0
        
    
    aberrationsurvival = np.exp(-aber)
    #print(aber)
    base_survival = np.exp(-(prebreak*defaultSurvParams['baseRate']))
    
    ret = apoptoticsurvival*aberrationsurvival*base_survival
    #print(ret, type(ret))
    return ret
        
    
    
    
    
    
    
def g2_survival(aber, postbreak):
    '''
							number of DSB at mitosis: post repair
							mitotic rate
	'''
	
	
	
	
    mitoticbreak = min(postbreak, 20)
	
    mitoticsurvival = np.exp(-(mitoticbreak * defaultSurvParams['mitoticRate']))
    
    aberrationsurvival = np.exp(-aber)
    
    ret = mitoticsurvival*aberrationsurvival
    #print(ret, type(ret))
    return ret
    
    
    
    
    
    
    
def m_survival(aber, prebreak, postbreak):
    
    mitoticbreak = postbreak
    
    mitoticsurvival = np.exp(-(mitoticbreak * defaultSurvParams['mitoticRate']))
    
    g1survival = g1_survival(aber, prebreak, postbreak)
    
    ret = mitoticsurvival*g1survival
    #print(ret, type(ret))
    return ret
    
    
    
    
    


def get_survival_rates(g1arrest):
    
    LET = LET_value
    dose = dose_value
    
    all_phase_rate = []
    
    for t in range(len(cell_phase)):
        all_phase_rate.append(calculatesurvival(cell_phase[t]))
        
        #if LET == 0 or particle == 0:
        #    all_phase_rate.append(calculatesurvival(cell_phase[t]))
            #print(all_phase_rate)
        #else:
        #    all_phase_rate.append(survivalmod(cell_phase[t], LET, dose))
        
    return all_phase_rate



def process_data():
    '''
    phase_survival is average of all runs in a file, separated in 3 phases [g1, g2, m]
    should theoretically be accurate for our NICE model since its composed of only 1 run per SDD
    
    async_survival is average of all 3 phases, list length = amount of runs in sdd
    '''
    sr = get_survival_rates(functional_g1arrest)
    #sys.stdout.write(str(sr))
    #sys.stdout.flush()
    all_rate = np.array(sr, np.int32)
    #print(all_rate)
    
    if plot_survival_rate == True:
        
        for i in all_rate:
            
            hist, bin_edges = np.histogram(i, bins=np.logspace(np.log10(1E-10), np.log10(1.0), 11) )
            bin_center = (bin_edges[:-1] + bin_edges[1:])/2
            
            plt.xscale("log")
            plt.plot(bin_center, hist)
            plt.xlabel("survival rate")
            plt.ylabel("count")
            plt.show()
        
    phase_survival = []
    for x in all_rate:
        mea = np.mean(x)
        phase_survival.append(mea)
        
    #phase_survival = [np.mean(x) for x in all_rate]
    #print("The phase survival rate is: ", phase_survival)
    
    new_all_rate = all_rate.transpose() #3 by 100, all_rate[0] is rate for all 3 phases for the first index
    
    async_rate = [np.mean(x) for x in new_all_rate]
    async_rate = np.mean(async_rate)
    
    #print("The overall survival rate is: ", async_rate)

    return phase_survival, async_rate


'''def survivalmod(phase, LET, dose):
    
    #if LET == 0 or particle == 0:
    #    return survival condition based on damage only
    

    
    
    DSBPerMBPPerGy = 35.0/6100.0
    DSBPerGray = DSBPerMBPPerGy*6100
    DSBPerkeV = DSBPerGray/(6.2415*4.0/3.0*np.pi*pow(defaultDNAParams['rNuc'],3))
    DSBPerUm = DSBPerkeV * LET
    rNuc_cell = pow(defaultCell["dna"]/6100.0, 1.0/3.0) * defaultDNAParams['rNuc']
    DSBPerTrack = (4.0/3.0) * rNuc_cell * DSBPerUm
    dosePerTrack = DSBPerTrack / (defaultCell["dna"]*DSBPerMBPPerGy)
    tracksPerNucleus = dose / dosePerTrack
    
    
    if DSBPerTrack<0.5: return calculatesurvival(phase)
    
    
    minTracks = max(1,int(tracksPerNucleus-3*np.sqrt(tracksPerNucleus)-3))
    maxTracks = int(tracksPerNucleus+3*np.sqrt(tracksPerNucleus)+3)
    trackRange = list(range(minTracks,maxTracks+1))
    trackWeights = poisson.pmf(trackRange,tracksPerNucleus)
    
    
    for n,weight in zip(trackRange,trackWeights):
        singlesurvival = np.array(calculatesurvival(phase)) * weight
        #print(n)
        tpn =[np.exp(-tracksPerNucleus)]*len(singlesurvival)
                
        survTotal = tpn + singlesurvival
    
    return survTotal'''


def statistics(repeats):
    '''
    takes in an integer > 0 to run the simulation for "repeats" amount of times
    '''

    phase_rate = []
    async_rate = []
    rate_and_std = []
    
    for run in range(repeats):
        
        #sys.stdout = sys.__stdout__
        #print(run)
        #sys.stdout.flush()
        
        phase_r, asyc_r = process_data()
        
        phase_rate.append(phase_r)
        async_rate.append(asyc_r)
        
        
    #print(phase_rate)
    phase_rate = np.transpose(phase_rate)
    

    for x in phase_rate:
        if plot_phase_rate == True:
            plt.hist(x)
            plt.show()
        
        rate_and_std.append([np.mean(x), np.std(x)])
        #print(np.mean(x), np.std(x))
        
    if plot_phase_rate == True:    
        plt.hist(async_rate) #async rate need to be recoded for contribution?
        plt.show()
        
    
    return rate_and_std








result = statistics(repeats)
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



xax = np.linspace(0, 4)
plt.yscale("log")
plt.plot(xax, linq1(xax), label = "chaudhary14, AG01522 cells")
plt.plot(xax, linq2(xax), label = "chaudhary14, U-87 cells")
plt.errorbar(1, result[0][0],yerr = result[0][1], fmt='r.', label = "this work")
plt.legend()
plt.title("Model Survival Rate In Comparison To PIDE Data")
plt.xlabel("Dose (Gy)")
plt.ylabel("Survival Rate")
plt.show()
