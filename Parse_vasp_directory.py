import os
import json
from chgnet.utils import parse_vasp_dir

# Base directory containing K/Li/Na → Materials → Structures
base_dir = 'path to base directory'
final_dataset_dicts = []

# Loop over first-level directories
for first_level in os.listdir(base_dir):
    first_level_path = os.path.join(base_dir, first_level)
    if os.path.isdir(first_level_path) and not first_level.startswith('.'):
        # Loop over second-level directories
        for second_level in os.listdir(first_level_path):
            second_level_path = os.path.join(first_level_path, second_level)
            if os.path.isdir(second_level_path) and not second_level.startswith('.'):
                # Loop over third-level directories
                for third_level in os.listdir(second_level_path):
                    third_level_path = os.path.join(second_level_path, third_level)
                    if os.path.isdir(third_level_path) and not third_level.startswith('.'):
                        try:
                            print(f"Parsing VASP data from: {third_level_path}")
                            dataset_dict = parse_vasp_dir(
                                base_dir=third_level_path,
                                save_path=os.path.join(third_level_path, "chgnet_dataset_5.json")
                            )

                            num_steps = len(dataset_dict["structure"])

                            # Every 5th step only
                            for i in range(0, num_steps, 5):
                                final_data = {
                                    'structure': dataset_dict['structure'][i].as_dict(),
                                    'uncorrected_total_energy': dataset_dict['uncorrected_total_energy'][i],
                                    'energy_per_atom': dataset_dict['energy_per_atom'][i],
                                    'force': dataset_dict['force'][i],
                                    'stress': dataset_dict['stress'][i] if dataset_dict.get('stress') else None,
                                    'magmom': dataset_dict['magmom'][i] if dataset_dict.get('magmom') else None,
                                }
                                final_dataset_dicts.append(final_data)

                            print(f"Added {num_steps // 5 + 1} data points from {third_level_path}")

                        except RuntimeError as e:
                            print(f"Error parsing {third_level_path}: {e}")

# Save the aggregated dataset
final_dataset_path = os.path.join(base_dir, "final_chgnet_dataset_5.json")
with open(final_dataset_path, "w") as f:
    json.dump(final_dataset_dicts, f)

print(f"Aggregated dataset saved to {final_dataset_path}")

