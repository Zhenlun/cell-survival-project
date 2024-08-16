import matplotlib.pyplot as plt
%matplotlib inline
import re
import os
os.chdir(os.getcwd())
import numpy as np
import csv
import math
from scipy.optimize import curve_fit
from scipy.stats import linregress


# Function to read data from a file and return lists
def read_data(file_name):
    survival_rate = []
    standard_deviation = []
    real_dose = []
    target_dose = []
    energies = []
    dsb_count = []
    aberration_count = []
    aberration_std_dev = []
    aberration_std_error = []
    
    with open(file_name, 'r') as file:
        for line in file:
            values = eval(line)
            survival_rate.append(float(values[0]))
            standard_deviation.append(float(values[1]))
            real_dose.append(float(values[2]))
            
            # Extract the energy from the file path
            path_parts = os.path.normpath(values[3]).split(os.sep)
            energy_part = path_parts[-1] if len(path_parts[-1].split('_')) > 1 else path_parts[-2]
            energy = energy_part.split('_')[-1]  # get the last part after splitting by '_'
            tdose = energy_part.split('_')[0] 
            if energy.replace('.', '', 1).isdigit():
                energies.append(float(energy))  # convert to float and append
                target_dose.append(float(tdose))
            else:
                raise Exception('Energy not found in file path!')
            
            # Append the DSB and aberration counts
            dsb_count.append(float(values[4]))
            aberration_count.append(float(values[5]))
            aberration_std_dev.append(float(values[6]))
            aberration_std_error.append(float(values[7]))
    
    #print(target_dose)            
    return survival_rate, standard_deviation, real_dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, target_dose

def filter_data_by_energy(survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, energy_value):
    filtered_indices = [i for i, e in enumerate(energies) if e == energy_value]
    
    filtered_survival_rate = [survival_rate[i] for i in filtered_indices]
    filtered_standard_deviation = [standard_deviation[i] for i in filtered_indices]
    filtered_dose = [dose[i] for i in filtered_indices]
    filtered_energies = [energies[i] for i in filtered_indices]
    filtered_dsb_count = [dsb_count[i] for i in filtered_indices]
    filtered_aberration_count = [aberration_count[i] for i in filtered_indices]
    filtered_aberration_std_dev = [aberration_std_dev[i] for i in filtered_indices]
    filtered_aberration_std_error = [aberration_std_error[i] for i in filtered_indices]
    
    return filtered_survival_rate, filtered_standard_deviation, filtered_dose, filtered_energies, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_dev, filtered_aberration_std_error

def filter_data_by_target_dose(survival_rate, standard_deviation, dose, target_dose, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, target_dose_value):
    filtered_indices = [i for i, td in enumerate(target_dose) if td == target_dose_value]
    
    filtered_survival_rate = [survival_rate[i] for i in filtered_indices]
    filtered_standard_deviation = [standard_deviation[i] for i in filtered_indices]
    filtered_dose = [dose[i] for i in filtered_indices]
    filtered_target_dose = [target_dose[i] for i in filtered_indices]
    filtered_dsb_count = [dsb_count[i] for i in filtered_indices]
    filtered_aberration_count = [aberration_count[i] for i in filtered_indices]
    filtered_aberration_std_dev = [aberration_std_dev[i] for i in filtered_indices]
    filtered_aberration_std_error = [aberration_std_error[i] for i in filtered_indices]
    
    return filtered_survival_rate, filtered_standard_deviation, filtered_dose, filtered_target_dose, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_dev, filtered_aberration_std_error


def convert_to_standard_error(standard_deviation, N):
    standard_error = [sd / math.sqrt(N) for sd in standard_deviation]
    return standard_error

def sr(a, b, dose):
    return np.exp(-a * dose - b * dose**2)

def exp_func(x, a):
    return np.exp(-a * x)

def plot_combined_lines(data_list, dose, survival_rate, standard_error, dsb_count, aberration_count, aberration_std_error, selection, ax, target_dose, box=False):
    # Unpack the list
    energy, a_list, b_list, cellline = data_list

    # Generate x values for the lines
    x = np.linspace(min(dose), max(dose), 100)

    # Convert dose, dsb_count, aberration_count, and target_dose to numpy arrays
    dose = np.array(dose)
    dsb_count = np.array(dsb_count)
    aberration_count = np.array(aberration_count)
    target_dose = np.array(target_dose)

    # Plot the data points with error bars based on selection
    if 'survival_rate' in selection:
        ax.errorbar(dose, survival_rate, yerr=standard_error, fmt='o', label=f'{energy} MeV')
        ax.set_yscale("log")
        ax.set_ylim(0.01, 1)

    if 'dsb_count' in selection:
        if box:
            print("Creating box plots for DSB count")
            unique_doses = sorted(set(target_dose))
            print("Unique doses:", unique_doses)
            data = [dsb_count[np.where(target_dose == ud)] for ud in unique_doses]
            positions = [np.mean(dose[np.where(target_dose == ud)]) for ud in unique_doses]
            std_devs = [np.std(dose[np.where(target_dose == ud)]) for ud in unique_doses]
            print("Data for box plots:", data)
            print("Positions for box plots:", positions)
            print("Standard deviations of positions:", std_devs)
            ax.boxplot(data, positions=positions, widths=4*np.array(std_devs))
        else:
            ax.plot(dose, dsb_count, 's', label=f'{energy} MeV', alpha = 0.3)
            slope = np.linalg.lstsq(dose[:, np.newaxis], dsb_count, rcond=None)[0]
            y_fit = slope * dose
            ax.plot(dose, y_fit, '--', label=f'fit: slope = {slope[0]:.2f} DSB/Gy')

    if 'aberration_count' in selection:
        if box:
            print("Creating box plots for aberration count")
            unique_doses = sorted(set(target_dose))
            print("Unique doses:", unique_doses)
            data = [aberration_count[np.where(target_dose == ud)] for ud in unique_doses]
            positions = [np.mean(dose[np.where(target_dose == ud)]) for ud in unique_doses]
            std_devs = [np.std(dose[np.where(target_dose == ud)]) for ud in unique_doses]
            print("Data for box plots:", data)
            print("Positions for box plots:", positions)
            print("Standard deviations of positions:", std_devs)
            ax.boxplot(data, positions=positions, widths=4*np.array(std_devs))
        else:
            ax.errorbar(dose, aberration_count, yerr=aberration_std_error, fmt='v', label=f'{energy} MeV', alpha = 0.3)
            print(np.linalg.lstsq(dose[:, np.newaxis], aberration_count, rcond=None))
            slope = np.linalg.lstsq(dose[:, np.newaxis], aberration_count, rcond=None)[0]
            y_fit = slope * dose
            ax.plot(dose, y_fit, '--', label=f'fit: slope = {slope[0]:.2f} CA/Gy')

    # Plot a line for each cell line
    if 'survival_rate' in selection:
        popt, pcov = curve_fit(exp_func, dose, survival_rate, p0=(1))
        y_fit = exp_func(x, popt[0])
        ax.plot(x, y_fit, 'b--', label=f'fit: slope = {popt[0]:.2f}')
        
        for i in range(len(a_list)):
            y = sr(a_list[i], b_list[i], x)
            ax.plot(x, y, label=f'{cellline[i]} {energy} MeV')

def process_and_plot_combined_lines(file_name, data_list, N, selection, ax, box=False):
    # Read data from file
    survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, target_dose = read_data(file_name)
    
    # Filter data by energy
    energy = data_list[0]
    filtered_survival_rate, filtered_standard_deviation, filtered_dose, filtered_energies, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_dev, filtered_aberration_std_error = filter_data_by_energy(
        survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, energy)
    
    # Filter the corresponding target_dose values
    filtered_target_dose = [target_dose[i] for i in range(len(target_dose)) if energies[i] == energy]
    
    # Convert standard deviation to standard error
    filtered_standard_error = convert_to_standard_error(filtered_standard_deviation, N)
    
    # Initialize lists to store positions and standard deviations
    all_positions = []
    all_std_devs = []
    
    # Plot the data
    if box:
        unique_target_doses = sorted(set(filtered_target_dose))
        for td in unique_target_doses:
            filtered_td_survival_rate, filtered_td_standard_deviation, filtered_td_dose, filtered_td_target_dose, filtered_td_dsb_count, filtered_td_aberration_count, filtered_td_aberration_std_dev, filtered_td_aberration_std_error = filter_data_by_target_dose(
                filtered_survival_rate, filtered_standard_deviation, filtered_dose, filtered_target_dose, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_dev, filtered_aberration_std_error, td)
            
            plot_combined_lines(data_list, filtered_td_dose, filtered_td_survival_rate, filtered_standard_error, filtered_td_dsb_count, filtered_td_aberration_count, filtered_td_aberration_std_error, selection, ax, filtered_td_target_dose, box)
            
            # Store positions and standard deviations
            positions = [np.mean(filtered_td_dose)]
            std_devs = [np.std(filtered_td_dose)]
            all_positions.extend(positions)
            all_std_devs.extend(std_devs)
            # Set x-ticks and labels
            ax.set_xticks(all_positions)
            ax.set_xticklabels([f'{pos:.2f}' for pos in all_positions])
        print("real dose mean:", all_positions)
        print("Standard deviations of dose mean:", all_std_devs)
    else:
        plot_combined_lines(data_list, filtered_dose, filtered_survival_rate, filtered_standard_error, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_error, selection, ax, filtered_target_dose, box)
    



particle = 'proton'
box = True
#1.38mev, tumor HeLa
e = 1.38
a = [0.53, 0.67]
b = [0.084, 0.037]
cellline = ["HeLa", "HeLaS3"]
list138 = [e, a, b, cellline]

#0.88mev, order is lung normal async HF-19, tongue tumor async SCC-25, peripheral blood normal G0 Lymphocytes
e = 0.88
a = [0.52, 0.81, 1.34]
b = [-0.053, -0.023, -0.467]
cellline = ["HF-19", "SCC-25", "G0 lymphocytes"]
list088 = [e, a, b, cellline]

#5.04mev, order iso lung normal async HF-19, tongue tumor async SCC-25
e = 5.04
a = [0.55, 0.41]
b = [0.074, 0.092]
cellline = ["HF-19", "SCC-25"]
list504 = [e, a, b, cellline]

# Create two figures
fig1, ax1 = plt.subplots(figsize=(6, 5))
fig2, ax2 = plt.subplots(figsize=(6, 5))

#data_file_name = "cluster_multiple_10000.txt"
data_file_name = "cluster_single_10000.txt"
N = 10000

# First plot with three lines
selection = ["aberration_count"]
process_and_plot_combined_lines(data_file_name, list088, N, selection, ax1, box=False)
process_and_plot_combined_lines(data_file_name, list088, N, selection, ax1, box=box)

#process_and_plot_combined_lines(data_file_name, list138, N, selection, ax1, box=box)
#process_and_plot_combined_lines(data_file_name, list504, N, selection, ax1, box=box)

# Last plot with three lines
selection = ["dsb_count"]
process_and_plot_combined_lines(data_file_name, list088, N, selection, ax2, box=False)
#process_and_plot_combined_lines(data_file_name, list088, N, selection, ax2, box=box)

#process_and_plot_combined_lines(data_file_name, list138, N, selection, ax2, box=box)
#process_and_plot_combined_lines(data_file_name, list504, N, selection, ax2, box=box)

# Add labels and titles
ax1.set_xlabel('Dose (Gy)', fontsize = 16)
ax1.set_ylabel('CA Count', fontsize = 16)
ax1.set_title('Single Cell Chromosome Aberration Count', fontsize = 16)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(True)
ax1.legend()
ax1.set_xlim([0, None])
ax1.set_ylim([0, None])
ax1.grid()

ax2.set_xlabel('Dose (Gy)', fontsize = 16)
ax2.set_ylabel('DSB Count', fontsize = 16)
ax2.set_title('Single Cell DSB Count', fontsize = 16)
ax2.legend()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(True)
ax2.spines['left'].set_visible(True)
ax2.legend()
ax2.set_xlim([0, None])
ax2.set_ylim([0, None])
ax2.grid()
# Show the plots
plt.show()

