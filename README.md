# Automated-electrode-screening
This repository contains the automated pipeline for the determination of voltage profile for 2D MT2 type material with the help of fine tuned CHGNet model for Li, Na, and K ion batteries. The only input is unit cell of 2D MT2 type material and the output will be voltage profile for Li, Na, and K. The scheme of the automated pipeline is as follows.

![image](https://github.com/user-attachments/assets/18c79082-80a8-4981-9a9e-6fef06a4571a)


# Prerquisites
Installation of pymatgen and CHGNet

The fine tune of pre-trained CHGNet has done following the CHGNet github page. The link is below
https://github.com/CederGroupHub/chgnet/tree/main/examples

All the materials has been extracted from the aNANt database. The link is below

https://anant.mrc.iisc.ac.in/apps/2D


