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
    
    print(target_dose)            
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


def convert_to_standard_error(standard_deviation, N):
    standard_error = [sd / math.sqrt(N) for sd in standard_deviation]
    return standard_error

def sr(a, b, dose):
    return np.exp(-a * dose - b * dose**2)

def exp_func(x, a):
    return np.exp(-a * x)
'''
def plot_data(data_list, dose, survival_rate, standard_error, dsb_count, aberration_count, aberration_std_error, selection):
    # Unpack the list
    energy, a_list, b_list, cellline = data_list

    # Generate x values for the lines
    x = np.linspace(min(dose), max(dose), 100)

    # Create a figure and a set of subplots
    fig, ax1 = plt.subplots()

    # Convert dose, dsb_count and aberration_count to numpy arrays
    dose = np.array(dose)
    dsb_count = np.array(dsb_count)
    aberration_count = np.array(aberration_count)

    # Plot the data points with error bars based on selection
    if 'survival_rate' in selection:
        ax1.errorbar(dose, survival_rate, yerr=standard_error, fmt='o', label='Survival Rate')
        ax1.set_ylabel('Survival Rate')
        ax1.set_yscale("log")
        ax1.set_ylim(0.01, 1)
        #ax1.legend(loc='upper left')

    if 'dsb_count' in selection:
        if ax1.get_ylabel():
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            ax2.plot(dose, dsb_count, 's', label='DSB Count')
            ax2.set_ylabel('DSB Count')  # we already handled the x-label with ax1
            #ax2.legend(loc='upper right')

            # Perform linear regression
            slope = np.linalg.lstsq(dose[:, np.newaxis], dsb_count, rcond=None)[0]
            y_fit = slope * dose
            plt.plot(dose, y_fit, '--', label=f'fit: y={slope[0]:.3f}x')
        else:
            ax1.plot(dose, dsb_count, 's', label='DSB Count')
            ax1.set_ylabel('DSB Count')
            slope = np.linalg.lstsq(dose[:, np.newaxis], dsb_count, rcond=None)[0]
            y_fit = slope * dose
            ax1.plot(dose, y_fit, '--', label=f'fit: y={slope[0]:.3f}x')
            #ax1.legend(loc='upper left')

    if 'aberration_count' in selection:
        if ax1.get_ylabel():
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            ax2.errorbar(dose, aberration_count, yerr=aberration_std_error, fmt='v', label='Aberration Count')
            ax2.set_ylabel('Aberration Count')  # we already handled the x-label with ax1
            #ax2.legend(loc='upper right')

            # Perform linear regression
            slope = np.linalg.lstsq(dose[:, np.newaxis], aberration_count, rcond=None)[0]
            y_fit = slope * dose
            ax2.plot(dose, y_fit, '--', label=f'fit: y={slope[0]:.3f}x')
        else:
            ax1.errorbar(dose, aberration_count, yerr=aberration_std_error, fmt='v', label='lethal aberrations')
            ax1.set_ylabel('Aberration Count')
            slope = np.linalg.lstsq(dose[:, np.newaxis], aberration_count, rcond=None)[0]
            y_fit = slope * dose
            ax1.plot(dose, y_fit, '--', label=f'fit: y={slope[0]:.3f}x')
            #ax1.legend(loc='upper left')
    


    # Plot a line for each cell line
    if 'survival_rate' in selection:
        popt, pcov = curve_fit(exp_func, dose, survival_rate, p0=(1))
        y_fit = exp_func(x, popt[0])
        ax1.plot(x, y_fit, 'b--', label='fit: a=%5.3f' % tuple(popt))
        
        for i in range(len(a_list)):
            y = sr(a_list[i], b_list[i], x)
            ax1.plot(x, y, label=cellline[i])
    
    # Add labels and title
    ax1.set_xlabel('Dose')
    plt.title(f'{particle} {energy} MeV')

    ax1.legend()
    if len(selection) > 1:
        ax2.legend()

    # Show the plot
    plt.show()


def process_and_plot_data(file_name, data_list, N, selection):
    # Read data from file
    survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, target_dose = read_data(file_name)
    
    # Filter data by energy
    energy = data_list[0]
    filtered_survival_rate, filtered_standard_deviation, filtered_dose, filtered_energies, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_dev, filtered_aberration_std_error = filter_data_by_energy(survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, energy)
    
    # Convert standard deviation to standard error
    filtered_standard_error = convert_to_standard_error(filtered_standard_deviation, N)
    
    # Plot the data
    plot_data(data_list, filtered_dose, filtered_survival_rate, filtered_standard_error, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_error, selection)
'''    
    
def plot_combined_lines(data_list, dose, survival_rate, standard_error, dsb_count, aberration_count, aberration_std_error, selection, ax):
    # Unpack the list
    energy, a_list, b_list, cellline, plotparticle = data_list
        
    # Generate x values for the lines
    x = np.linspace(min(dose), max(dose), 100)

    # Convert dose, dsb_count and aberration_count to numpy arrays
    dose = np.array(dose)
    dsb_count = np.array(dsb_count)
    aberration_count = np.array(aberration_count)

    # Plot the data points with error bars based on selection
    if 'survival_rate' in selection:
        ax.errorbar(dose, survival_rate, yerr=standard_error, fmt='o', label=f'{energy} MeV {plotparticle}')
        ax.set_yscale("log")
        ax.set_ylim(0.01, 1)

    if 'dsb_count' in selection:
        #ax.errorbar(dose, dsb_count, fmt='s', label=f'{energy} MeV {plotparticle}')
        ax.plot(dose, dsb_count, 's', label=f'{energy} MeV {plotparticle}')
        slope = np.linalg.lstsq(dose[:, np.newaxis], dsb_count, rcond=None)[0][0]
        residuals = dsb_count - (slope * dose)
        SE_slope = np.sqrt(np.sum(residuals**2) / (len(dose) - 1)) / np.sqrt(np.sum(dose**2))
        y_fit = slope * dose
        ax.plot(dose, y_fit, '--', label=f'fit: slope = {slope:.2f} ± {SE_slope:.2f} DSB/Gy')

    if 'aberration_count' in selection:
        #ax.errorbar(dose, aberration_count, yerr=aberration_std_error, fmt='v', label=f'{energy} MeV {plotparticle}')
        ax.plot(dose, aberration_count, 'v', label=f'{energy} MeV {plotparticle}')
        slope = np.linalg.lstsq(dose[:, np.newaxis], aberration_count, rcond=None)[0][0]
        residuals = aberration_count - (slope * dose)
        SE_slope = np.sqrt(np.sum(residuals**2) / (len(dose) - 1)) / np.sqrt(np.sum(dose**2))
        y_fit = slope * dose
        ax.plot(dose, y_fit, '--', label=f'fit: slope = {slope:.3f} ± {SE_slope:.3f} CA/Gy')

    # Plot a line for each cell line
    if 'survival_rate' in selection:
        popt, pcov = curve_fit(exp_func, dose, survival_rate, p0=(1))
        y_fit = exp_func(x, popt[0])
        ax.plot(x, y_fit, 'b--', label=f'fit: slope = {popt[0]:.2f}')
        
        for i in range(len(a_list)):
            y = sr(a_list[i], b_list[i], x)
            ax.plot(x, y, label=f'{cellline[i]} {energy} MeV')

def process_and_plot_combined_lines(file_name, data_list, N, selection, ax):
    # Read data from file
    survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, target_dose = read_data(file_name)
    
    # Filter data by energy
    energy = data_list[0]
    filtered_survival_rate, filtered_standard_deviation, filtered_dose, filtered_energies, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_dev, filtered_aberration_std_error = filter_data_by_energy(survival_rate, standard_deviation, dose, energies, dsb_count, aberration_count, aberration_std_dev, aberration_std_error, energy)
    
    # Convert standard deviation to standard error
    filtered_standard_error = convert_to_standard_error(filtered_standard_deviation, N)
    
    # Plot the data
    plot_combined_lines(data_list, filtered_dose, filtered_survival_rate, filtered_standard_error, filtered_dsb_count, filtered_aberration_count, filtered_aberration_std_error, selection, ax)

    
"-------------------------------------------------------------------"

particle = 'proton'
box = True
#1.38mev, tumor HeLa
e = 1.38
a = [0.53, 0.67]
b = [0.084, 0.037]
cellline = ["HeLa", "HeLaS3"]
list138 = [e, a, b, cellline, particle]

#0.88mev, order is lung normal async HF-19, tongue tumor async SCC-25, peripheral blood normal G0 Lymphocytes
e = 0.88
a = [0.52, 0.81, 1.34]
b = [-0.053, -0.023, -0.467]
cellline = ["HF-19", "SCC-25", "G0 lymphocytes"]
list088 = [e, a, b, cellline, particle]

#5.04mev, order iso lung normal async HF-19, tongue tumor async SCC-25
e = 5.04
a = [0.55, 0.41]
b = [0.074, 0.092]
cellline = ["HF-19", "SCC-25"]
list504 = [e, a, b, cellline, particle]

"-------------------------------------------------------------------"


'''selection = ("survival_rate", "dsb_count")
data_file_name = "cluster_single_100.txt"
N = 100
process_and_plot_data(data_file_name, list088, N, selection)
process_and_plot_data(data_file_name, list138, N, selection)
process_and_plot_data(data_file_name, list504, N, selection)'''

'''
data_file_name = "cluster_multiple_10000.txt"
N = 10000
#selection = ("survival_rate", "aberration_count")
selection = ["aberration_count"]
process_and_plot_data(data_file_name, list088, N, selection)
process_and_plot_data(data_file_name, list138, N, selection)
process_and_plot_data(data_file_name, list504, N, selection)


#selection = ("survival_rate", "dsb_count")
selection = ["dsb_count"]
process_and_plot_data(data_file_name, list088, N, selection)
process_and_plot_data(data_file_name, list138, N, selection)
process_and_plot_data(data_file_name, list504, N, selection)'''

"-------------------------------------------------------------------"

# Create two figures
fig1, ax1 = plt.subplots(figsize=(6, 5))
fig2, ax2 = plt.subplots(figsize=(6, 5))

data_file_name = "cluster_multiple_10000.txt"
#data_file_name = "single_10000.txt"
N = 10000

# First plot with three lines
selection = ["aberration_count"]
process_and_plot_combined_lines(data_file_name, list088, N, selection, ax1)
process_and_plot_combined_lines(data_file_name, list138, N, selection, ax1)
process_and_plot_combined_lines(data_file_name, list504, N, selection, ax1)

# Last plot with three lines
selection = ["dsb_count"]
process_and_plot_combined_lines(data_file_name, list088, N, selection, ax2)
process_and_plot_combined_lines(data_file_name, list138, N, selection, ax2)
process_and_plot_combined_lines(data_file_name, list504, N, selection, ax2)

data_file_name = "felixphoton.txt"
particle = 'X-ray'
list025 = [0.25, a, b, cellline, particle]
selection = ["aberration_count"]
process_and_plot_combined_lines(data_file_name, list025, N, selection, ax1)
selection = ["dsb_count"]
process_and_plot_combined_lines(data_file_name, list025, N, selection, ax2)



#now for reference
dose = np.array([0, 0.25, 0.5, 0.75, 1, 2, 3, 3.5])
aberration_count = np.array([1, 0.0383, 0.0595, 0.0980, 0.1388, 0.3540, 0.6239, 1.0575])

ax1.errorbar(dose, aberration_count, yerr=[0, 0.0125, 0.0174, 0.0240, 0.0245, 0.0410, 0.0519, 0.1089], 
             fmt='r.', label="in-vitro X-ray")

# Calculate the slope
slope = np.linalg.lstsq(dose[:, np.newaxis], aberration_count, rcond=None)[0][0]
residuals = aberration_count - (slope * dose)
SE_slope = np.sqrt(np.sum(residuals**2) / (len(dose) - 1)) / np.sqrt(np.sum(dose**2))
y_fit = slope * dose

ax1.plot(dose, y_fit, 'r', label=f'fit: slope = {slope:.2f} ± {SE_slope:.2f} CA/Gy')



# Add labels and titles
ax1.set_xlabel('Dose (Gy)', fontsize = 16)
ax1.set_ylabel('CA Count', fontsize = 16)
ax1.set_title('100 cells combined Chromosome Aberration count', fontsize = 16)
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
ax2.set_title('100 cells combined DSB Count', fontsize = 16)
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

