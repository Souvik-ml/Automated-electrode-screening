import os

# Define a mapping of elements to their specific potential file
element_potcar_map = {
    'Ba': 'Ba_sv',
    'Be': 'Be',
    'Ca': 'Ca_pv',
    'Cd': 'Cd',
    'Co': 'Co',
    'Cu': 'Cu',
    'Fe': 'Fe',
    'Mg': 'Mg',
    'Mn': 'Mn',
    'Ni': 'Ni',
    'Pb': 'Pb',
    'Pd': 'Pd',
    'Pt': 'Pt',
    'Sr': 'Sr_sv',
    'Zn': 'Zn',
    'Ce': 'Ce',
    'Ge': 'Ge',
    'Hf': 'Hf',
    'Ir': 'Ir',
    'Mo': 'Mo',
    'Os': 'Os',
    'Re': 'Re',
    'Ru': 'Ru',
    'Sn': 'Sn',
    'Tc': 'Tc',
    'Te': 'Te',
    'Ti': 'Ti',
    'W': 'W',      
    'Zr': 'Zr_sv_GW',
    'Br': 'Br',
    'C': 'C',
    'N': 'N',
    'O': 'O',
    'Cl': 'Cl',
    'F': 'F',
    'I': 'I',
    'S': 'S',
    'H': 'H',
    'P': 'P',
    'Se': 'Se',
    'Li': 'Li',
    'Na': 'Na',
    'K': 'K_pv'  
}

# Define the base directory where POTCAR files are stored
potential_base_dir = 'path to POTENTIAL files'

# Define the main path where the structure directories are located
main_path = 'path to structures'

# Loop through directories and subdirectories to find POSCAR files
for folder in os.listdir(main_path):
    folder_path = os.path.join(main_path, folder)
    if os.path.isdir(folder_path):  # Check if it's a directory
        for sub_folder in os.listdir(folder_path):
            sub_folder_path = os.path.join(folder_path, sub_folder)
            if os.path.isdir(sub_folder_path):  # Check if it's a directory
                for sub_folder_2 in os.listdir(sub_folder_path):
                    sub_folder_2_path = os.path.join(sub_folder_path, sub_folder_2)
                    if os.path.isdir(sub_folder_2_path):  # Check if it's a directory
                        for filename in os.listdir(sub_folder_2_path):
                            if filename == 'POSCAR':  # If the file is POSCAR
                                poscar_path = os.path.join(sub_folder_2_path, filename)

                                # Read the POSCAR file to get the list of elements
                                with open(poscar_path, 'r') as f:
                                    lines = f.readlines()

                                # Line 6 of POSCAR contains the element symbols
                                elements = lines[5].split()

                                # Prepare the POTCAR file by concatenating individual element POTCARs
                                output_potcar_path = os.path.join(sub_folder_2_path, 'POTCAR')  # Path for the combined POTCAR file

                                # Open the final POTCAR file for writing
                                with open(output_potcar_path, 'wb') as potcar:
                                    for element in elements:
                                        # Remove any extra spaces/newlines and ensure element is a string
                                        element = element.strip()

                                        # Get the specific potential file for this element
                                        if element not in element_potcar_map:
                                            raise KeyError(f"No potential mapping found for element: {element}")

                                        specific_potential = element_potcar_map[element]

                                        # Construct the path to the desired POTCAR file
                                        element_potcar_path = os.path.join(potential_base_dir, specific_potential, 'POTCAR')

                                        # Check if the POTCAR file exists
                                        if not os.path.exists(element_potcar_path):
                                            raise FileNotFoundError(f"POTCAR file for element {element} not found at {element_potcar_path}")

                                        # Read and append the POTCAR file content
                                        with open(element_potcar_path, 'rb') as elem_potcar:
                                            potcar.write(elem_potcar.read())

                                print(f"POTCAR file successfully created at: {output_potcar_path}")

