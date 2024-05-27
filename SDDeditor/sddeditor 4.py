import os
import re

def update_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.startswith('Data entries,'):
                parts = line.split(',')
                if parts[7].strip() == '0':  # Check if the 7th number is 0
                    parts[7] = ' 1'  # Replace with 1 (note the space for formatting)
                line = ','.join(parts)
            file.write(line)

def search_and_update_files(start_path):
    for root, dirs, files in os.walk(start_path):
        for file in files:
            if file.endswith('SDDOutput.txt'):
                update_file(os.path.join(root, file))

# Use the directory of the script file as the start directory
start_directory = os.path.dirname(os.path.realpath(__file__))
search_and_update_files(start_directory)
