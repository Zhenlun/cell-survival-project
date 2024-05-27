import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Replace 'your_data.csv' with the path to your CSV file
data = pd.read_csv('PIDE.csv')

# Function to calculate survival rate
def survival_rate(a, b, dose):
    return np.exp(-a * dose - b * dose * dose)

# Plotting function with ion type filter
def plot_lq_model(data, variable, ion_type, variable_range=None, rows=None, start=None, end=None):
    # Filter data for the specified ion type
    data = data[data['Ion'] == ion_type]
    
    # Filter data for the specified variable range if provided
    if variable_range is not None:
        min_val, max_val = variable_range
        data = data[(data[variable] >= min_val) & (data[variable] <= max_val)]
    
    if rows is not None:
        data = data[data['#ExpID'].isin(rows)]
    elif start is not None and end is not None:
        data = data[(data['#ExpID'] >= start) & (data['#ExpID'] <= end)]
    
    # Create a figure and a set of subplots
    fig, ax = plt.subplots()
    
    # Generate a colormap index based on LET or Energy values
    cmap = plt.cm.viridis
    norm = plt.Normalize(data[variable].min(), data[variable].max())
    
    for i, row in data.iterrows():
        a = row['ai_paper'] if pd.notna(row['ai_paper']) else row['ai_fit']
        b = row['bi_paper'] if pd.notna(row['bi_paper']) else row['bi_fit']
        doses = np.linspace(0, 5, 100)  # Adjust dose range as needed
        survival_rates = survival_rate(a, b, doses)
        
        # Map the LET or Energy value to a color
        color = cmap(norm(row[variable]))
        
        #ax.plot(doses, survival_rates, label=f'Exp {row["#ExpID"]}', color=color)
        ax.plot(doses, survival_rates, color=color, alpha = 0.1)
    
    # Set the y-axis to logarithmic scale
    ax.set_yscale('log')
    
    # Add a color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label(variable)
    
    ax.set_xlabel('Dose')
    ax.set_ylabel('Survival Rate')
    ax.set_title(f'Linear Quadratic Model of Cell Survival for {ion_type} Ions')
    ax.legend()
    plt.show()

# Example usage:
plot_lq_model(data, 'LET', '1H')
#plot_lq_model(data, 'Energy', '1H')
plot_lq_model(data, 'DNAcontent', '1H')
#plot_lq_model(data, 'LET', '1H', variable_range=(0, 10))
plot_lq_model(data, 'Energy', '1H', variable_range=(0, 70))
#plot_lq_model(data, 'DNAcontent', '1H', variable_range=(0, 10))
#plot_lq_model(data, 'LET', '4He', variable_range=(10, 20))
# plot_lq_model(data, 'LET', '4He', rows=[1, 2, 3])  # For specific rows and ion type '4He'
# plot_lq_model(data, 'LET', '12C', start=1, end=100)  # For a range of rows and ion type '12C'