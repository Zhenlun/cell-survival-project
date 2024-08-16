import os
import shutil
import re

def move_files(base_dir):
    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)
        if os.path.isdir(subdir_path):
            for filename in os.listdir(subdir_path):
                filepath = os.path.join(subdir_path, filename)
                if os.path.isfile(filepath):
                    # Split the filename into name and extension
                    name, ext = os.path.splitext(filename)
                    # Split the name by underscore
                    parts = name.split('_')
                    # The index is the last part before the extension
                    index = parts[-1]
                    index = ''.join(filter(str.isdigit, index))
                    # Check if the index is a number
                    if index != "" and index.isdigit():
                        new_dir = os.path.join(subdir_path, index)
                        os.makedirs(new_dir, exist_ok=True)
                        shutil.move(filepath, new_dir)
                        print(f'Moved file {filename} from folder {subdir} to folder {new_dir}')

base_dir = os.getcwd()
move_files(base_dir)
