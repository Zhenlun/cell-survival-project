import os
import shutil

def process_data(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # Read the header section
        line = f_in.readline()
        while line.strip() != "***EndOfHeader***;" and line:
            f_out.write(line)  # Write the header line to the output file
            
            if line.startswith("Chromosome sizes"):
                chromosome_sizes = list(map(float, line.strip(";\n").split(",")[2:]))
            line = f_in.readline()
        if line:  # Check if the end of header line was found
            f_out.write("***EndOfHeader***;\n")  # Write the end of header line to the output file
            # Process the data section
            for line in f_in:
                fields = line.split(";")
                index = int(fields[2].split(",")[1])  # Get the index from field 3
                chromosome_size = chromosome_sizes[index-1] * 1e6  # Convert to basepairs
                position = int(fields[3])  # Get the chromosome position from field 4
                new_position = position / chromosome_size  # Calculate the new ratio
                fields[3] = str(new_position)  # Replace the original position data
                f_out.write(";".join(fields))  # Write the modified data to the new file

# Walk through each subfolder and process all the data files
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith("_SDDOutput.txt"):
            path = os.path.join(root, file)
            temp_file = os.path.join(root, 'temp.txt')
            process_data(path, temp_file)
            shutil.move(temp_file, path)