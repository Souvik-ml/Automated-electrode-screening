from chgnet.model import CHGNet
import numpy as np
from pymatgen.core import Structure

from chgnet.utils import parse_vasp_dir, read_json

from chgnet.data.dataset import StructureData, get_train_val_test_loader


dataset_list = read_json("path to json file where the extracted structures have been stored")
from chgnet.data.dataset import StructureData, get_train_val_test_loader
structures = [Structure.from_dict(entry["structure"]) for entry in dataset_list]
energies_per_atom = [entry["energy_per_atom"] for entry in dataset_list]
forces = [entry["force"] for entry in dataset_list]
stresses = [entry.get("stress") for entry in dataset_list]
magmoms = [entry.get("magmom") for entry in dataset_list]

dataset = StructureData(
    structures=structures,
    energies=energies_per_atom,
    forces=forces,
    stresses=stresses,  # can be None
    magmoms=magmoms,  # can be None
)
train_loader, val_loader, test_loader = get_train_val_test_loader(
    dataset, batch_size=256, train_ratio=0.9, val_ratio=0.05
)

from chgnet.model import CHGNet
from chgnet.trainer import Trainer


# Load pretrained CHGNet
chgnet = CHGNet.load()



# Optionally fix the weights of some layers
for layer in [
    chgnet.atom_embedding,
    chgnet.bond_embedding,
    chgnet.angle_embedding,
    chgnet.bond_basis_expansion,
    chgnet.angle_basis_expansion,
    chgnet.atom_conv_layers[:-1],
    chgnet.bond_conv_layers,
    chgnet.angle_layers,
]:
    for param in layer.parameters():
        param.requires_grad = False

# Define Trainer
trainer = Trainer(
    model=chgnet,
    targets="efs",
    optimizer="Adam",
    scheduler="CosLR",
    criterion="MSE",
    epochs=500,
    learning_rate=1e-4,
    use_device="cpu",
    print_freq=6,
)

trainer.train(train_loader, val_loader, test_loader)

model = trainer.model
best_model = trainer.best_model  # best model based on validation energy MAE




