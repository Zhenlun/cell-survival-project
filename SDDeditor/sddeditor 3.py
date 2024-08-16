import os
import re

def sum_numbers_in_file(filename):
    with open(filename, 'r') as file:
        numbers = file.read().split()
        return sum(float(num) for num in numbers if num.replace('.', '', 1).isdigit())
    
def update_file(file_path, real_dose):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.startswith('Dose or fluence,'):
                parts = line.split(',')
                parts[-1] = str(real_dose)+';\n'
                if len(parts) == 2 and parts[1].strip().replace(';', '').replace('.', '', 1).isdigit():
                    line = f'Dose or fluence, 1, {parts[1]}'
            file.write(line)
    print(f"Modified file: {file_path}")

def search_and_update_files(start_path):
    for root, dirs, files in os.walk(start_path):
        for file in files:
            if file.startswith('SDDOutput'):
                # Find the corresponding "RealDoseDep" file in the same directory
                for real_dose_file in files:
                    if real_dose_file.startswith('RealDoseDep'):
                        real_dose = sum_numbers_in_file(os.path.join(root, real_dose_file))
                        break  # We found the file, no need to check the rest
                else:
                    continue
                
                update_file(os.path.join(root, file), real_dose)

# Use the directory of the script file as the start directory
start_directory = os.path.dirname(os.path.realpath(__file__))
search_and_update_files(start_directory)
