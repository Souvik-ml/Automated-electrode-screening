import os
import shutil

# Define the base directory where your subdirectories are located
BASE_DIR = "/home/souvik/barc/database"

# Define the target directory where filtered structures will be copied
DEST_DIR = "/home/souvik/barc/homo_genious_2D_mat"
os.makedirs(DEST_DIR, exist_ok=True)

# Iterate through all items in the base directory
for subdir in os.listdir(BASE_DIR):
    subdir_path = os.path.join(BASE_DIR, subdir)
    
    # Check if the subdirectory is a directory
    if os.path.isdir(subdir_path):
        
        # Split the directory name into components (e.g., 'NO-Cl-Cd-2D' -> ['NO', 'Cl', 'Cd', '2D'])
        parts = subdir.split('-')
        
        # Check if the first two parts (atom1 and atom2) are the same
        if parts[0] == parts[1]:  # atom1 == atom2
            print(f"Matching directory: {subdir_path}")
            
            # Look for the 'poscar' directory inside the matching directory
            poscar_dir = os.path.join(subdir_path, "poscars")
            
            if os.path.isdir(poscar_dir):
                print(f"Found poscar directory: {poscar_dir}")
                
                # Copy all structure files in 'poscar' to the destination
                for file_name in os.listdir(poscar_dir):
                    file_path = os.path.join(poscar_dir, file_name)
                    if os.path.isfile(file_path):
                        # Debug: Print file being copied
                        print(f"Copying file: {file_path}")
                        
                        shutil.copy(file_path, os.path.join(DEST_DIR, file_name))
            else:
                print(f"'poscar' directory not found in: {subdir_path}")
        else:
            print(f"Skipping directory (atom1 != atom2): {subdir_path}")

print(f"Filtered structures have been copied to {DEST_DIR}.")



