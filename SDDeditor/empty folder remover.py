import os

# Set the directory you want to start from
root_folder = os.getcwd()

# Walk through all subdirectories and files in the directory
for folder_path, subfolders, files in os.walk(root_folder, topdown=False):
    # Check if the directory is empty
    if not subfolders and not files:
        try:
            os.rmdir(folder_path)
            print(f'Deleted empty folder: {folder_path}')
        except OSError as e:
            print(f'Error: {folder_path} : {e.strerror}')
