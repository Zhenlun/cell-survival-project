import os

# Set the directory you want to start from
root_folder = os.getcwd()

# Walk through all subdirectories and files in the directory
for folder_path, subfolders, files in os.walk(root_folder):
    # Loop through each file in the current directory
    for filename in files:
        # Check if the file is a .txt file
        if filename.endswith('.txt'):
            # Construct full file path
            file_path = os.path.join(folder_path, filename)
            # Check if the file is empty
            if os.stat(file_path).st_size == 0:
                # Remove the file
                os.remove(file_path)
                print(f'Deleted empty file: {filename}')
