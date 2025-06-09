from pymatgen.core import Structure, Element
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from chgnet.model import StructOptimizer
import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import random
from chgnet.model import CHGNet

input_folder_path = 'path to unit cell of 2D material structures'

current_path = os.getcwd()

chgnet_ft = CHGNet.from_file("path to the best fine tuned model or tar file generated during the training of CHGNet model")

adsorbed_structures_path = os.path.join(current_path, 'Adsorbed_structures')
os.makedirs(adsorbed_structures_path, exist_ok=True)

optimized_structure_path = os.path.join(current_path, 'Optimized_structures')
os.makedirs(optimized_structure_path, exist_ok=True)

trajectory_path = os.path.join(current_path, 'Trajectory_path')
os.makedirs(trajectory_path, exist_ok=True)

ion_dirs = ['Li', 'Na', 'K']


for ion in ion_dirs:
  os.makedirs(os.path.join(adsorbed_structures_path, ion), exist_ok=True)
  os.makedirs(os.path.join(optimized_structure_path, ion), exist_ok=True)
  os.makedirs(os.path.join(trajectory_path, ion), exist_ok=True)


for structure in os.listdir(input_folder_path):
  for ion in ion_dirs:
    working_ion_path = os.path.join(adsorbed_structures_path, ion)
    adsorb_path = os.path.join(working_ion_path, structure)
    os.makedirs(adsorb_path, exist_ok=True)

    # Load the structure
    struct = Structure.from_file(os.path.join(input_folder_path, structure))

    # Activate the asf function
    asf = AdsorbateSiteFinder(struct)
    sites = asf.find_adsorption_sites()
    adsorption_sites = {key:value for key, value in sites.items() if key!='all'}
    # Iterate through the adsorbate sites
    for site_name, site_coords in adsorption_sites.items():
      for idx, site_coord in enumerate(site_coords):
        # Create a new structure by adding the adsorbate at the given site
        new_opt_struct = struct.copy()
        # Convert site_coord from Cartesian to fractional coordinates
        site_coord = new_opt_struct.lattice.get_fractional_coords(site_coord)
        # Append the Li atom at the current site
        new_opt_struct.append(Element(ion), site_coord)
        # Define file names for the current position (index-based naming)
        file_name_cif = os.path.join(adsorb_path, f"{structure}_{site_name}_{idx}.cif")
        #file_name_poscar = f"POSCAR_{site_name}_{idx+1}"
        # Save the structure in CIF format
        new_opt_struct.to(fmt="cif", filename=file_name_cif)
        # Save the structure in POSCAR format
        #new_opt_struct.to(fmt="poscar", filename=file_name_poscar)
        print(f"Saved structure with {ion} at {site_name} position {idx} to {file_name_cif}")
        #print(f"Saved structure with Li at {site_name} position {idx} to {file_name_cif} and {file_name_poscar}")


for ion in ion_dirs:
  for folder in os.listdir(adsorbed_structures_path):
    if ion in folder:
      folder_path = os.path.join(adsorbed_structures_path, folder)

      for sub_folder in os.listdir(folder_path):

        sub_folder_path = os.path.join(folder_path, sub_folder)

        adsorb_ion_struct_path = os.path.join(optimized_structure_path, ion)
        opt_struct_path = os.path.join(adsorb_ion_struct_path, sub_folder)
        os.makedirs(opt_struct_path, exist_ok=True)

        # Create subdirectories for saving trajectory paths
        adsorb_opt_ion_traj_path = os.path.join(trajectory_path, ion)
        adsorb_opt_traj_path = os.path.join(adsorb_opt_ion_traj_path, sub_folder)
        os.makedirs(adsorb_opt_traj_path, exist_ok=True)


        for filename in os.listdir(sub_folder_path):

          struc_path = os.path.join(sub_folder_path, filename)

          base_name = os.path.splitext(filename)[0]



          opt = StructOptimizer(model = chgnet_ft)
          opt_result = opt.relax(Structure.from_file(struc_path), fmax=0.01, save_path = os.path.join(adsorb_opt_traj_path, f'{base_name}.pickle'))
          opt_result['final_structure'].to(fmt='cif', filename=os.path.join(opt_struct_path, filename))


most_stable_adsorb_structures_path = 'path to most_stable_adsorb_structures'
os.makedirs(most_stable_adsorb_structures_path, exist_ok=True)

op_struct_dict = {}

for dir in os.listdir(trajectory_path):
  dir_path = os.path.join(trajectory_path, dir)
  ion_dir_path = os.path.join(most_stable_adsorb_structures_path, dir)
  os.makedirs(ion_dir_path, exist_ok=True)

  for sub_dir in os.listdir(dir_path):
    sub_dir_path = os.path.join(dir_path, sub_dir)

    ion_sub_dir_path = os.path.join(ion_dir_path, sub_dir)
    os.makedirs(ion_sub_dir_path, exist_ok=True)

    Energy = {}

    for files in os.listdir(sub_dir_path):
      if files.endswith('.pickle'):
        with open(os.path.join(sub_dir_path, files), 'rb') as f:
          result = pickle.load(f)
          Energy[files] = result['energy'][-1]

        # Save the Energy dictionary in the subdir as a pickle file (optional)
    energy_pickle_path = os.path.join(sub_dir_path, 'Energy.pickle')
    with open(energy_pickle_path, 'wb') as f:
      pickle.dump(Energy, f)


    # Find the file with the minimum energy
    min_energy_file = min(Energy, key=Energy.get)
    min_energy_value = Energy[min_energy_file]

    # Prepare to copy the corresponding CIF file
    min_energy_file_base_name = f'{os.path.splitext(min_energy_file)[0]}.cif'
    source_cif_path = os.path.join(optimized_structure_path, dir, sub_dir, min_energy_file_base_name)
    target_cif_path = os.path.join(ion_sub_dir_path, min_energy_file_base_name)

    if os.path.exists(source_cif_path):
      shutil.copy(source_cif_path, target_cif_path)
      print(f"Copied {min_energy_file_base_name} to {ion_sub_dir_path}")
    else:
      print(f"File {min_energy_file_base_name} not found in {os.path.join(optimized_structure_path, dir, sub_dir)}")

    # Update the op_struct_dict with the minimum energy file and value
    op_struct_dict[min_energy_file] = min_energy_value
    print(f"The minimum energy file is {min_energy_file} with energy {min_energy_value}")

# Save the op_struct_dict as a pickle file in the current directory
op_struct_dict_path = 'path to op_struct_dict.pickle'
with open(op_struct_dict_path, 'wb') as f:
  pickle.dump(op_struct_dict, f)

print(f"'op_struct_dict' has been saved as a pickle file at: {op_struct_dict_path}")


# Create the folder for supercell structures
supercell_path = 'path to most_stable_adsorb_structures_supercell'
os.makedirs(supercell_path, exist_ok=True)

# Loop through all files in most_stable_Li_adsorb_structures
for dir in os.listdir(most_stable_adsorb_structures_path):
    # Create a subdirectory for each structure
    input_dir_path = os.path.join(most_stable_adsorb_structures_path, dir)

    dir_path = os.path.join(supercell_path, dir)
    os.makedirs(dir_path, exist_ok=True)

    for sub_dir in os.listdir(input_dir_path):

      input_subdir_path = os.path.join(input_dir_path, sub_dir)
      subdir_path = os.path.join(dir_path, sub_dir)
      os.makedirs(subdir_path, exist_ok=True)

      for filename in os.listdir(input_subdir_path):
        base_name = os.path.splitext(filename)[0]
        # Load the structure and create a supercell
        structure = Structure.from_file(os.path.join(input_subdir_path, filename))
        supercell_structure = structure.copy()
        supercell_structure.make_supercell([2, 2, 1])

        # Start with the full Li structure
        mod_filename = f'{base_name}_4{dir}.cif'
        supercell_structure.to(fmt='cif', filename=os.path.join(subdir_path, mod_filename))

        # Iteratively remove Li atoms to create reduced-Li structures
        for n in range(3, -1, -1):  # 3Li, 2Li, 1Li, 0Li
          n_ion = int(supercell_structure.composition[dir])

          if n_ion > 0:  # Remove a Li atom if there are any left
            remove_ids = random.sample(list(range(n_ion)), 1)
            supercell_structure.remove_sites(remove_ids)

          # Save the structure for each step, including when Li is 0
          mod_filename = f'{base_name}_{n}{dir}.cif'
          supercell_structure.to(fmt='cif', filename=os.path.join(subdir_path, mod_filename))

print(f"Supercell structures saved in {supercell_path}")


# Optimize the supercell structures with 0Li, 1Li, 2Li, 3Li, and 4Li
optimized_supercell_path = 'path to Optimized_supercell'
os.makedirs(optimized_supercell_path, exist_ok=True)

opt_super_traj_path = 'path to Optimized_supercell_trajectory'
os.makedirs(opt_super_traj_path, exist_ok=True)

for dir in os.listdir(supercell_path):

  input_dir_path = os.path.join(supercell_path, dir)

  # Create subdirectory to save the optimized structures within Optimized_supercell directory
  opt_dir_path = os.path.join(optimized_supercell_path, dir)
  os.makedirs(opt_dir_path, exist_ok=True)

  # Create directories and subdirectories to save the trajectory files
  opt_dir_traj_path = os.path.join(opt_super_traj_path, dir)
  os.makedirs(opt_dir_traj_path, exist_ok=True)

  for sub_dir in os.listdir(input_dir_path):
    input_sub_dir_path = os.path.join(input_dir_path, sub_dir)

    opt_sub_dir_path = os.path.join(opt_dir_path, sub_dir)
    os.makedirs(opt_sub_dir_path, exist_ok=True)

    opt_sub_dir_traj_path = os.path.join(opt_dir_traj_path, sub_dir)
    os.makedirs(opt_sub_dir_traj_path, exist_ok=True)

    for structure in os.listdir(input_sub_dir_path):
      struc_path = os.path.join(input_sub_dir_path, structure)
      base_name = os.path.splitext(structure)[0]
      struc = Structure.from_file(struc_path)
      opt = StructOptimizer(model=chgnet_ft)
      opt_result = opt.relax(struc, fmax=0.01, save_path = os.path.join(opt_sub_dir_traj_path, f'{base_name}.pickle'))

      # Save optimized structure
      optimized_structure_path = os.path.join(opt_sub_dir_path, f'{base_name}.cif')
      opt_result['final_structure'].to(fmt='cif', filename=optimized_structure_path)

      print(f"Optimized structure of {structure} saved to: {optimized_structure_path}")


# Track the Energy of each step of lithiation for each Material
for dir in os.listdir(opt_super_traj_path):

  input_dir_path = os.path.join(opt_super_traj_path, dir)

  for sub_dir in os.listdir(input_dir_path):
    input_sub_dir_path = os.path.join(input_dir_path, sub_dir)

    Energy = {}

    # Process all pickle files in the subdir
    for files in os.listdir(input_sub_dir_path):
      files_base_name = os.path.splitext(files)[0]
      if files.endswith('.pickle'):
        with open(os.path.join(input_sub_dir_path, files), 'rb') as f:
          result = pickle.load(f)
          Energy[files_base_name] = result['energy'][-1]

    # Save the Energy dictionary in the subdir as a pickle file (optional)
    energy_pickle_path = os.path.join(input_sub_dir_path, 'Energy.pickle')
    with open(energy_pickle_path, 'wb') as f:
      pickle.dump(Energy, f)


# Directory to save voltage dictionaries
voltage_dic_save_path = 'path to voltage_dic_save_folder'
os.makedirs(voltage_dic_save_path, exist_ok=True)


reference_energies = {
    "Li": -1.88,
    "Na": -1.29,
    'K': -1.09
}

# Iterate through each element (Li, Na, K)
for element, reference_energy in reference_energies.items():
  element_voltage_path = os.path.join(voltage_dic_save_path, element)
  os.makedirs(element_voltage_path, exist_ok=True)

  # Loop through all directories in the trajectory path
  for dir in os.listdir(opt_super_traj_path):
    if element in dir:
      input_dir_path = os.path.join(opt_super_traj_path, dir)

      for sub_dir in os.listdir(input_dir_path):
        input_sub_dir_path = os.path.join(input_dir_path, sub_dir)

        for filename in os.listdir(input_sub_dir_path):
          if filename.startswith('Energy'):
            # Load the energy dictionary
            with open(os.path.join(input_sub_dir_path, filename), 'rb') as f:
              voltage_energy = pickle.load(f)

              # Sort the energy dictionary based on the number of ions
              sorted_voltage_energy = dict(sorted(voltage_energy.items(), key=lambda item: int(item[0].split('_')[-1].replace(dir, ''))))

              E_Voltage_list = list(sorted_voltage_energy.values())
              voltage_dict = {}

              # Calculate voltages
              for number_of_ions in range(0, 4):
                voltage_dict[f'V{number_of_ions}{number_of_ions + 1}'] = -(E_Voltage_list[number_of_ions + 1] - (E_Voltage_list[number_of_ions] + reference_energy))

        # Save the voltage dictionary
        with open(os.path.join(element_voltage_path, f'{sub_dir}.pickle'), 'wb') as f:
          pickle.dump(voltage_dict, f)

        print(f"Voltage dictionary saved for {sub_dir} of {element} at {os.path.join(element_voltage_path, f'{sub_dir}.pickle')}")


### Combine plot

combine_voltage_plot_save_path = 'path to voltage_plot_save_folder_combined'
os.makedirs(combine_voltage_plot_save_path, exist_ok=True)

# Get list of materials common to all folders
material_files = {ion: set(os.listdir(os.path.join(voltage_dic_save_path, ion))) for ion in ion_dirs}
common_files = set.intersection(*material_files.values())  # Find common files for all ions

# Process each common material
for filename in common_files:
    if filename.endswith('.pickle'):
        try:
            # Load voltage data for Li and Na
            voltage_data = {}
            for ion in ion_dirs:
                file_path = os.path.join(voltage_dic_save_path, ion, filename)
                with open(file_path, 'rb') as f:
                    voltage_data[ion] = pickle.load(f)

            # Extract the material descriptor (e.g., "O-O-Ti")
            descriptor = filename.split('_')[1].replace('.pickle', '')

            # Split the descriptor into atoms (e.g., ['O', 'O', 'Ti'])
            split_name = descriptor.split('-')

            # Assign atoms correctly
            atom1 = split_name[2]  # "Ti"
            atom2 = split_name[0]  # "O"

            # Generate the material name
            material_name = f"{atom1}({atom2})₂"

            # Prepare data for plotting
            number_of_ions = list(range(-1, 5))
            number_of_ions.remove(0)  # Remove 0
            ion_annotations = []
            plot_data = {}
            for ion in ion_dirs:
                voltage_list = list(voltage_data[ion].values())
                voltage_list.insert(0, voltage_list[0])  # Duplicate first value for the plot
                plot_data[ion] = voltage_list

                # Calculate Average Voltage and OCV
                avg_voltage = sum(voltage_data[ion].values()) / len(voltage_data[ion])
                ocv = voltage_list[-1]
                ion_annotations.append(f"{ion} ({avg_voltage:.2f}V, {ocv:.2f}V)")

            # Create combined plot
            plt.figure(figsize=(8, 6))
            for ion in ion_dirs:
                plt.step(number_of_ions, plot_data[ion], 'o-', linewidth=3, markersize=10, label=f'{ion} ion')

            min_voltage = min([min(plot_data[ion]) for ion in ion_dirs])
            max_voltage = max([max(plot_data[ion]) for ion in ion_dirs])

            # Annotate plot
            plt.xlabel('Number of Ions', fontsize=14, fontweight='bold')
            plt.ylabel('Voltage (V)', fontsize=14, fontweight='bold')
            plt.text(0.5, 0.9, material_name, fontsize=16, fontweight='bold', transform=plt.gca().transAxes, ha='center')

            plt.text(
                0.98, 0.95,
                "\n\n".join(ion_annotations),
                fontsize=12, fontweight='bold', transform=plt.gca().transAxes,
                ha='right', va='top', bbox=dict(facecolor='white', alpha=0.8, edgecolor='black')
            )
            plt.legend(fontsize=12, loc='upper left', frameon=True, prop={'weight': 'bold'})
            plt.xlim(0, 4.5)
            plt.ylim(min_voltage - 0.5, max_voltage + 0.7)
            plt.xticks(fontsize=12, fontweight='bold')
            plt.yticks(fontsize=12, fontweight='bold')

            # Save combined plot
            plot_path = os.path.join(combine_voltage_plot_save_path, f"{material_name}_Li_vs_Na_vs_K_voltage.png")
            plt.savefig(plot_path, dpi=500, bbox_inches='tight')
            plt.clf()
            #print(f"Combined plot saved: {plot_path}")

        except Exception as e:
            print(f"Error processing material {filename}: {e}")


