
import os
import math
import numpy as np
from scipy import stats
import io
import sys

from concurrent.futures import ThreadPoolExecutor

from repairanalysis import medrasrepair
from repairanalysis import sddparser
from threading import Lock



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

def run_simulation(sddp):
    # Create an in-memory text stream
    memory_file = io.StringIO()
    
    # Redirect stdout to the in-memory stream
    original_stdout = sys.stdout
    sys.stdout = memory_file

    # Run the simulation
    medrasrepair.repairSimulation(sddp, 'Spectrum')

    # Reset stdout
    sys.stdout = original_stdout
    
    # Reset the stream position to the beginning
    memory_file.seek(0)
    
    # Read the data from the in-memory stream
    repair_data = [line.split() for line in memory_file if line[0].isdigit()]
    
    # Close the in-memory stream
    memory_file.close()
    
    return repair_data

def record_viability(num_simulations, sddpath, single):
    if not single:
        data_cell_count = run_simulation(sddpath)
        simulations_to_run = math.ceil(num_simulations / len(data_cell_count))
    else:
        simulations_to_run = num_simulations
    viability_records = []
    dsb_records = []
    aber_records = []
        
    for _ in range(simulations_to_run):
        repair_data = run_simulation(sddpath)
        for data_set in repair_data:
            # Check if the last element is "misrepair!", indicating a viability of 1
            if data_set[-1] == "misrepair!":
                viability_records.append('1')
                dsb_records.append(data_set[1])
                aber_records.append("0")
            else:
                # Otherwise, record the actual last element as the viability
                viability_records.append(data_set[-1])
                dsb_records.append(data_set[1])
                aber_records.append(data_set[-2])

    return np.array(viability_records, dtype=int), np.array(dsb_records, dtype=int), np.array(aber_records, dtype=int)
''' #multi_thread
lock = Lock()

def run_simulation(sddp):
    # Create an in-memory text stream
    memory_file = io.StringIO()
    
    with lock:
        # Redirect stdout to the in-memory stream
        original_stdout = sys.stdout
        sys.stdout = memory_file

        # Run the simulation
        medrasrepair.repairSimulation(sddp, 'Spectrum')

        # Reset stdout
        sys.stdout = original_stdout
    
    # Reset the stream position to the beginning
    memory_file.seek(0)
    
    # Read the data from the in-memory stream
    repair_data = [line.split() for line in memory_file if line[0].isdigit()]
    
    # Close the in-memory stream
    memory_file.close()
    
    return repair_data

def record_viability(num_simulations, sddpath, single):
    threads = 16
    if not single:
        data_cell_count = run_simulation(sddpath)
        simulations_to_run = math.ceil(num_simulations / len(data_cell_count))
    else:
        simulations_to_run = num_simulations
    viability_records = []
    dsb_records = []
    aber_records = []
        
    # Create a ThreadPoolExecutor with num_threads
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Use list comprehension to create a list of futures
        futures = [executor.submit(run_simulation, sddpath) for _ in range(simulations_to_run)]
        
        for future in futures:
            repair_data = future.result()
            for data_set in repair_data:
                # Check if the last element is "misrepair!", indicating a viability of 1
                if data_set[-1] == "misrepair!":
                    viability_records.append('1')
                    dsb_records.append(data_set[1])
                    aber_records.append("0")
                else:
                    # Otherwise, record the actual last element as the viability
                    viability_records.append(data_set[-1])
                    dsb_records.append(data_set[1])
                    aber_records.append(data_set[-2])

    return np.array(viability_records, dtype=int), np.array(dsb_records, dtype=int), np.array(aber_records, dtype=int)'''


def plot_histogram_and_calculate_survival_rate(viability_records):
    # Calculate the survival rate
    survival_rate = np.count_nonzero(viability_records == 1) / len(viability_records)
    return survival_rate

    
    

def calculate_survival_rate(num_simulations, sddpath, single):
    viability_results, dsb_array, aber_array = record_viability(num_simulations, sddpath, single)
    #print(len(viability_results))
    survival_rate = plot_histogram_and_calculate_survival_rate(viability_results)
    return survival_rate, dsb_array, aber_array

def calculate_statistics(number_of_trials, sddpath, single):
    survival_rate, dsb_array, aber_array = calculate_survival_rate(number_of_trials, sddpath, single)
    n = len(aber_array)
    p = survival_rate
    variance = p * (1 - p) / n
    standard_deviation = np.sqrt(variance)
    
    dsb = np.mean(dsb_array)
    dsb_std = np.std(dsb_array)
    dsb_se = stats.sem(dsb_array)
    
    aber = np.mean(aber_array)
    aber_std_dev = np.std(aber_array)
    aber_se = stats.sem(aber_array)
    return survival_rate, variance, standard_deviation, dsb, aber, aber_std_dev, aber_se
        

if __name__ == "__main__":
    path = os.getcwd()
    sddpath = path
    folder_rate = []
    folder_variance = []
    folder_deviation = []
    folder_dsb = []
    folder_aber = []
    folder_aber_std_dev = []
    folder_aber_se = []
    
    paths = []
    repeats = 10000
    '''allfoldername = ['_inner_electron', '_inter_electron', '_outer_electron', 
                             '_inner_proton', '_inter_proton', '_outer_proton', 
                             '_inner_alpha', '_inter_alpha', '_outer_alpha']
    allfoldername = ['gy']'''
    #allfoldername =['proton']
    allfoldername = [str(i) for i in range(1, 101)]
    excludefoldername = ["proton"]
    #excludefoldername = ["gasdfg"]
    
    single = True
    write_file_name = "test_10000.txt"
    
    with open(write_file_name, "r") as input_file:
        for line in input_file:
            if line.strip():
                data = eval(line.strip())
                paths.append(data[3])

    with open(write_file_name, "a") as output_file:
        for root, subfolders, files in os.walk(sddpath):
            for foldername in subfolders:
                if any(sub in foldername for sub in allfoldername) and any(sub not in foldername for sub in excludefoldername):
                    folder_path = os.path.join(root, foldername)
                    if folder_path in paths:
                        continue
                    rate, var, std_dev, dsb, aber, aber_std_dev, aber_se = calculate_statistics(repeats, folder_path, single)
                    folder_rate.append(rate)
                    folder_variance.append(var)
                    folder_deviation.append(std_dev)
                    folder_dsb.append(dsb)
                    folder_aber.append(aber)
                    folder_aber_std_dev.append(aber_std_dev)
                    folder_aber_se.append(aber_se)
                    
                    folder_dose = []
                    for froot, fsubfolders, ffiles in os.walk(folder_path):
                        for file in ffiles:
                            if file.endswith('.txt'):
                                file_path = os.path.join(folder_path, file)
                                header, events = sddparser.parseSDDFile(file_path)
                                folder_dose.append(header['Dose or Fluence'][1])
                    dose = np.mean(folder_dose)
                    print([rate, std_dev, dose, folder_path, dsb, aber, aber_std_dev, aber_se], file=output_file, flush=True)
    
    '''with open(write_file_name, "a") as output_file:
        for root, subfolders, files in os.walk(sddpath):
            for foldername in subfolders:
                if any(sub in foldername for sub in allfoldername) and any(sub not in foldername for sub in excludefoldername):
                    folder_path = os.path.join(root, foldername)
                    rate, var, std_dev, dsb, aber, aber_std_dev, aber_se = calculate_statistics(repeats, folder_path, single)
                    folder_rate.append(rate)
                    folder_variance.append(var)
                    folder_deviation.append(std_dev)
                    folder_dsb.append(dsb)
                    folder_aber.append(aber)
                    folder_aber_std_dev.append(aber_std_dev)
                    folder_aber_se.append(aber_se)
                    
                    folder_dose = []
                    for froot, fsubfolders, ffiles in os.walk(folder_path):
                            for file in ffiles:
                                if file.endswith('.txt'):
                                    file_path = os.path.join(folder_path, file)
                                    header, events = sddparser.parseSDDFile(file_path)
                                    folder_dose.append(header['Dose or Fluence'][1])
                    dose = np.mean(folder_dose)
                    print([rate, std_dev, dose, folder_path, dsb, aber, aber_std_dev, aber_se], file=output_file, flush=True)'''