
#place in damage_results_merged folder to work


import os



path = os.getcwd()
#sddpath = path + "\\examplesdd.txt"

fileNames = os.listdir(path)
filePaths = [path+"\\"+f for f in fileNames]

for files in filePaths:
    if "." in files or "__pycache__" in files:
        continue
    
    os.chdir(files)
    inside_data_folder = os.listdir()
    
    all_sdd_file_names = []
    for file in inside_data_folder:
        if "SDD" in file:
            all_sdd_file_names.append(file)

    for sdd_file in all_sdd_file_names:
        
        with open(sdd_file) as f:
        	#enable following lines 1 line at a time
        	#may be able to just do newText.replace() to replace subsequent lines but i didnt do any test
            newText = ""

            #newText=f.read().replace(";\nIrradiation target", ";\nDose rate, 0.0;\nIrradiation target")
            #newText=f.read().replace("***EndOfHeader***;\n1", "***EndOfHeader***;\n2")
            #newText=f.read().replace("variable ", "")
            #newText=f.read().replace("variable, variable ", "")
            #newText=f.read().replace("Non specified ", "0")
            #print(newText)
            f.close()
        
        with open(sdd_file, "w") as f:
            if newText != "":
                f.write(newText)
                f.close()

    
    
    
    
    
#header, event = s.parseSDDFile(sdd_file)
    