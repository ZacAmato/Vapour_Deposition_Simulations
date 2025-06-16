#!/usr/bin/python

import math
import matplotlib.pyplot as plt
import numpy as np
import os
import random
from scipy.stats import maxwell

gro_input_file = open('C:/Users/Gromacs/XY_PBC_Runs/Co-deposition/New_deposition/surface/same name/quick_silic_surface_1H2O_relax_20K.gro', 'r')

osef = gro_input_file.readline()

line_natoms = gro_input_file.readline()
line_natoms_cut = line_natoms.split()
natoms = int(line_natoms_cut[0])

x_atom = np.zeros(natoms)
y_atom = np.zeros(natoms)
z_atom = np.zeros(natoms)
v_x = np.zeros(natoms)
v_y = np.zeros(natoms)
v_z = np.zeros(natoms)
atom_type = np.zeros(natoms)
types = [None]*natoms
x_box = 0.0
y_box = 0.0
z_box = 0.0
nsilica = 12500
nwater = 0
z_ref = 0.0
ntot = 0

for iatom in range(natoms):

    line_to_read = gro_input_file.readline()
    line_cut = line_to_read.split()

    if iatom < 9999:

        types[iatom] = line_cut[1]
        x_atom[iatom] = float(line_cut[3])
        y_atom[iatom] = float(line_cut[4])
        z_atom[iatom] = float(line_cut[5])

        v_x[iatom] = float(line_cut[6])
        v_y[iatom] = float(line_cut[7])
        v_z[iatom] = float(line_cut[8])

        if iatom == 0:

            x_ref = x_atom[iatom]
            y_ref = y_atom[iatom]

    else:

        types[iatom] = line_cut[1][0:2]
        x_atom[iatom] = float(line_cut[2])
        y_atom[iatom] = float(line_cut[3])
        z_atom[iatom] = float(line_cut[4])

        v_x[iatom] = float(line_cut[5])
        v_y[iatom] = float(line_cut[6])
        v_z[iatom] = float(line_cut[7])

d_lash = 10.0

for jatom in range(natoms-4):

    dist = np.sqrt((x_atom[jatom]-x_atom[natoms-4])**2+(y_atom[jatom]-y_atom[natoms-4])**2+(z_atom[jatom]-z_atom[natoms-4])**2)

    if dist < d_lash:

        d_lash = dist

if d_lash < 1.0:
    
    print("YES")
    
else:
        
    print("NO")

    #os.system('cp output_md.gro input_addition.gro')

