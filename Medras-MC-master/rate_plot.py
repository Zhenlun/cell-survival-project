import matplotlib.pyplot as plt
%matplotlib inline
import re
import os
os.chdir(os.getcwd())
import numpy as np
import csv

set_energy = 0.25 #MeV

# Function to convert energy to eV
def energy_to_ev(energy_str):
    # Handle the case where energy is in the format '1-1MeV' which means 1.1 million eV
    if '-' in energy_str:
        base, decimal = energy_str.split('-')
        if 'MeV' in decimal:
            return (float(base) + float(decimal.replace('MeV', '')) / 10) * 1e6
        elif 'keV' in decimal:
            return (float(base) + float(decimal.replace('keV', '')) / 10) * 1e3
    elif 'keV' in energy_str:
        return float(energy_str.replace('keV', '')) * 1e3
    elif 'MeV' in energy_str:
        return float(energy_str.replace('MeV', '')) * 1e6
    elif 'eV' in energy_str:
        return float(energy_str.replace('eV', ''))
    else:
        return float(energy_str)

# Function to read data from a file and return lists
def read_data(file_name):
    survival_rate = []
    standard_deviation = []
    dose = []
    data = []
    energies = []
    
    with open(file_name, 'r') as file:
        for line in file:
            values = eval(line)
            survival_rate.append(float(values[0]))
            standard_deviation.append(float(values[1]))
            dose.append(float(values[2]))
            data.append(values[3])
            
            # Extract the energy from the file path
            path_parts = os.path.normpath(values[3]).split(os.sep)
            energy_part = path_parts[-1] if len(path_parts[-1].split('_')) > 1 else path_parts[-2]
            energy = energy_part.split('_')[-1]  # get the last part after splitting by '_'
            if energy.replace('.', '', 1).isdigit():
                energies.append(float(energy))  # convert to float and append
            else:
                raise Exception('Energy not found in file path!')
                
    return survival_rate, standard_deviation, dose, data, energies


'''def read_data(file_name):
    survival_rate = []
    standard_deviation = []
    dose = []
    data = []
    energies = []
    
    with open(file_name, 'r') as file:
        for line in file:
            values = eval(line)
            survival_rate.append(float(values[0]))
            standard_deviation.append(float(values[1]))
            dose.append(float(values[2]))
            data.append(values[3])
            if set_energy == 0:
                energy_strn = re.search(r'n(\d+[\.-]?\d*[kM]?eV)', values[3])
                energy_strx = re.search(r'x(\d+[\.-]?\d*[kM]?eV)', values[3])
                if energy_strn:
                    energy_strn = energy_strn.group(1)
                    energies.append(energy_to_ev(energy_strn))
                elif energy_strx:
                    energy_strx = energy_strx.group(1)
                    energies.append(energy_to_ev(energy_strx))
                else:
                    raise Exception('nothing found!')
            else:
                energies = [set_energy*1e6]*4
    return survival_rate, standard_deviation, dose, data, energies'''

# Filter the data based on particle type and location
def filter_data(data, dose, survival_rate, standard_deviation, energies, particle_type=None, location=None):
    filtered_indices = []
    title_parts = []

    if particle_type:
        title_parts.append(particle_type)
    else:
        title_parts.append('all particles')

    if location:
        title_parts.append(location)
    else:
        title_parts.append('all scoring volumes')

    title = 'Survival Rate for ' + ' and '.join(title_parts)
            
    for index, d in enumerate(data):
        if (particle_type is None or particle_type in d) and (location is None or location in d):
            filtered_indices.append(index)
    # Use list comprehension to create filtered lists based on the indices
    filtered_dose = [dose[i] for i in filtered_indices]
    filtered_survival_rate = [survival_rate[i] for i in filtered_indices]
    filtered_standard_deviation = [standard_deviation[i] for i in filtered_indices]
    filtered_energies = [energies[i] for i in filtered_indices]
    print((particle_type is None or particle_type in d) , (location is None or location in d))

    return filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, title

# Function to plot data
def plot_data(dose, survival_rate, standard_deviation, energies, title, a, b, color_map=plt.cm.viridis, label_suffix=''):
    min_energy = min(energies)
    max_energy = max(energies)
    for idx in range(len(dose)):
        energy_ev = energies[idx]
        if set_energy == 0:
            color = color_map((energy_ev - min_energy) / (max_energy - min_energy))
            #plt.errorbar(dose[idx], survival_rate[idx], yerr=standard_deviation[idx]**2, fmt='o', color=color, label=f'{energy_ev} eV {label_suffix}')
            plt.errorbar(dose[idx], survival_rate[idx], yerr=standard_deviation[idx]**2, fmt='or', color=colory_ev)
        else: 
            #plt.errorbar(dose[idx], survival_rate[idx], yerr=standard_deviation[idx]**2, fmt='o', label=f'{energy_ev} eV {label_suffix}')
            plt.errorbar(dose[idx], survival_rate[idx], yerr=standard_deviation[idx]**2, fmt='or')
    if set_energy == 0:
        plt.colorbar(plt.cm.ScalarMappable(cmap=color_map), label=f'Energy (eV) {label_suffix}')
    fitx = np.linspace(0, max(dose), 100)
    for i in range(len(a)):
        if len(name) != 0:
            plt.plot(fitx, sr(a[i], b[i], fitx), 'b', label = name[i])
        else:
            plt.plot(fitx, sr(a[i], b[i], fitx), 'b')
    plt.legend()
    
# Main function to process and plot files
def process_and_plot(file_name, a, b, particle_type=None, location=None, color_map=plt.cm.viridis, label_suffix=''):
    survival_rate, standard_deviation, dose, data, energies = read_data(file_name)
    filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, plot_title = filter_data(data, dose, survival_rate, standard_deviation, energies, particle_type, location)
    plot_data(filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, plot_title, a, b, color_map, label_suffix)
    return plot_title

def sr(a, b, dose):
    return np.exp(-a * dose - b * dose**2)

def process_csv(file_name):
    with open(file_name, 'r') as f:
        reader = csv.DictReader(f)
        expid_list = []
        alpha_list = []
        beta_list = []
        for row in reader:
            if row['Î±'] != 'N/A' and row['Î²'] != 'N/A':
                expid_list.append(int(row['#ExpID']))
                alpha_list.append(float(row['Î±']))
                beta_list.append(float(row['Î²']))
    return expid_list, alpha_list, beta_list

plt.figure()
particle = 'proton'
#a = [0.53, 0.67]
#b = [0.084, 0.037]
#1.38mev, tumor HeLa
a = [0.53]
b = [0.84]

#0.88mev, order is lung normal async HF-19, tongue tumor async SCC-25, peripheral blood normal G0 Lymphocytes
a = [0.52, 0.81, 1.34]
b = [-0.053, -0.023, -0.467]

#5.04mev, order iso lung normal async HF-19, tongue tumor async SCC-25
a = [0.55, 0.41]
b = [0l.074, 0.092]

#name = ['goodhead92, exp1', 'goodhead92, exp5']
#particle = 'photon'
#expid, alpha, beta = process_csv('RadPhysBio 250keV photon monoenergetic.csv')
#indexpre = 10
#indexpost = 20

#a = alpha[indexpre:indexpost]
#b = beta[indexpre:indexpost]
#name = expid[indexpre:indexpost]
#a = alpha
#b = beta
name = []


# Call the process_and_plot function for both new and old data
#title_new = process_and_plot(os.path.join(os.getcwd(), 'new_survival.txt'), particle_type=particle, color_map=plt.cm.Reds_r, label_suffix='(New)')
#title_old = process_and_plot(os.path.join(os.getcwd(), 'old_survival.txt'), particle_type=particle, color_map=plt.cm.Blues_r, label_suffix='(Old)')
#title_new = process_and_plot(os.path.join(os.getcwd(), 'new_survival.txt'), particle_type=particle, color_map=plt.cm.Reds, label_suffix='(New)')
#title_old = process_and_plot(os.path.join(os.getcwd(), 'old_survival.txt'), particle_type=particle, color_map=plt.cm.Blues, label_suffix='(Old)')
#title_new = process_and_plot(os.path.join(os.getcwd(), 'x0_25.txt'), a = a, b = b, particle_type=particle, color_map=plt.cm.Reds_r, label_suffix='(New)')
#itle_new = process_and_plot(os.path.join(os.getcwd(), 'p1_38.txt'), a = a, b = b, particle_type=particle, color_map=plt.cm.Reds_r, label_suffix='(New)')
title_new = process_and_plot(os.path.join(os.getcwd(), 'cluster_multiple.txt'), a = a, b = b, particle_type=particle, color_map=plt.cm.Reds_r)

# Finalize the plot

plt.ylim(0.001, 1.1)
plt.xlabel('Dose')
plt.ylabel('Survival Rate')
plt.title(f'{title_new}')
#plt.legend()
plt.show()

