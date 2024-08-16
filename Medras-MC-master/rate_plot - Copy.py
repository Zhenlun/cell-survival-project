import matplotlib.pyplot as plt
import re
import os

# Initialize lists for the data
survival_rate = []
standard_deviation = []
dose = []
data = []
energies = []  # List to store energy values in eV

file_namenew = os.path.join(os.getcwd(), 'new_survival.txt')
file_nameold = os.path.join(os.getcwd(), 'old_survival.txt')

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

# Read the data from the file
with open(file_name, 'r') as file:
    # Extract energy string and convert to eV
    for line in file:
        values = eval(line)
        survival_rate.append(float(values[0]))
        standard_deviation.append(float(values[1]))
        dose.append(float(values[2]))
        data.append(values[3])
        # Extract energy string and convert to eV
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


# Define min_energy and max_energy based on the extracted energies
min_energy = min(energies)
max_energy = max(energies)

# Filter the data based on particle type and location
def filter_data(particle_type=None, location=None):
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

    title = ' and '.join(title_parts)

            
    for index, d in enumerate(data):
        if (particle_type is None or particle_type in d) and (location is None or location in d):
            filtered_indices.append(index)
    # Use list comprehension to create filtered lists based on the indices
    filtered_dose = [dose[i] for i in filtered_indices]
    filtered_survival_rate = [survival_rate[i] for i in filtered_indices]
    filtered_standard_deviation = [standard_deviation[i] for i in filtered_indices]
    filtered_energies = [energies[i] for i in filtered_indices]

    title = ' and '.join(title_parts)
    return filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, title



# Example usage of filter_data function

filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, plot_title = filter_data(particle_type='proton')
#filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, plot_title = filter_data(particle_type='electron')
#filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, plot_title = filter_data(particle_type='alpha')
#filtered_dose, filtered_survival_rate, filtered_standard_deviation, filtered_energies, plot_title = filter_data()




# Color the data points based on energy
for idx in range(len(filtered_dose)):
    energy_ev = filtered_energies[idx]
    color = plt.cm.viridis((energy_ev - min_energy) / (max_energy - min_energy))
    #plt.errorbar(filtered_dose[idx], filtered_survival_rate[idx], yerr=filtered_standard_deviation[idx], fmt='o', color=color, label=f'{energy_ev} eV', ecolor='lightgray', elinewidth=3, capsize=0)
    plt.errorbar(filtered_dose[idx], filtered_survival_rate[idx], yerr=filtered_standard_deviation[idx], fmt='o', color=color, label=f'{energy_ev} eV')

# Set the colorbar to show the energy scale
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min_energy, vmax=max_energy))
plt.colorbar(sm, label='Energy (eV)')

plt.ylim(0.001, 1.1)
#plt.yscale('log')
plt.xlabel('Dose')
plt.ylabel('Survival Rate')
plt.title(f'Survival Rate for {plot_title}')
#plt.legend()
plt.show()

