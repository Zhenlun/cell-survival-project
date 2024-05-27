import os
import re

def update_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    with open(file_path, 'w') as file:
        for line in lines:
            if line.startswith('Dose or fluence,'):
                parts = line.split(',')
                if len(parts) == 2 and parts[1].strip().replace(';', '').replace('.', '', 1).isdigit():
                    line = f'Dose or fluence, 1,{parts[1]}'
            file.write(line)

def search_and_update_files(start_path):
    for root, dirs, files in os.walk(start_path):
        for file in files:
            if file.endswith('SDDOutput.txt'):
                update_file(os.path.join(root, file))

# Use the directory of the script file as the start directory
start_directory = os.path.dirname(os.path.realpath(__file__))
search_and_update_files(start_directory)
