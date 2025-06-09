import os
import pickle
import numpy as np

# Directory containing the pickle files
directory = 'path to metal ion present in voltage_dic_save_folder'

# List to store names of files that meet the condition
valid_files = []

# Iterate through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.pickle'):
        file_path = os.path.join(directory, filename)
        
        # Open and load the pickle file
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        
        # Convert dictionary values to a list
        value_ls = list(data.values())
        
        # Compute absolute differences between consecutive values
        diffs = [np.abs(value_ls[i] - value_ls[i+1]) for i in range(len(value_ls)-1)]
        
        # Check if all differences are less than 0.26
        if all(diff < 0.26 for diff in diffs):
            print(filename)  # Print file name
            valid_files.append(filename)  # Save file name

# Save the list of valid file names to a file
with open('valid_pickle_files_Li.txt', 'w') as f:
    for name in valid_files:
        f.write(name + '\n')

