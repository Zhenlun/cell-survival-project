import os

def modify_files(directory):
    for foldername, subfolders, filenames in os.walk(directory):
        for filename in filenames:
            if filename.startswith("SDDOutput"):
                file_path = os.path.join(foldername, filename)
                with open(file_path, 'r') as f:
                    content = f.read()
                content = content.replace(";\nIrradiation target", ";\nDose rate, 0.0;\nIrradiation target")
                content = content.replace("***EndOfHeader***;\n1", "***EndOfHeader***;\n2")
                content = content.replace("variable ", "")
                content = content.replace("variable, variable ", "")
                content = content.replace("Non specified", "0")
                content = content.replace("Time 0;", "Time, 0;")
                with open(file_path, 'w') as f:
                    f.write(content)
                print(f"Modified file: {file_path}")

modify_files(os.getcwd())
