#!/usr/bin/python

import math
import matplotlib.pyplot as plt
import numpy as np
import os
import random
from scipy.stats import maxwell

gro_file = str('300Kwater_on_20Ksilic_CO_addition')

output = str('300Kwater_on_20Ksilic_CO_addition')


gro_Zmax_file = open(f'../{gro_file}.gro', 'r')

truc = gro_Zmax_file.readline()
osef = gro_Zmax_file.readline()

Zmax = 0

for iline in range(200, 300):
    
    line = gro_Zmax_file.readline()
    line_cut = line.split()
        
    if float(line_cut[5]) > Zmax:
            
        Zmax = float(line_cut[3])
        
        

gro_input_file = open(f'../{gro_file}.gro', 'r')

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
nsilica = 0
nwater = 0
nco = 0
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

    if types[iatom] == 'SZ1':

        ntot = ntot+1
        nsilica = nsilica+1

        if z_ref < z_atom[iatom]:
            z_ref = z_atom[iatom]
            
    elif types[iatom] == 'SZ':
        
        ntot = ntot+1
        nsilica = nsilica+1
        
        if z_ref < z_atom[iatom]:
            z_ref = z_atom[iatom]

    elif types[iatom] == 'OW':

        if z_atom[iatom] - z_ref < 1.5 and z_atom[iatom] - Zmax < 3.0:
            ntot = ntot + 1
            nwater = nwater +1

            if z_ref < z_atom[iatom]:
                z_ref = z_atom[iatom]
                ntot = ntot + 1
                nwater = nwater +1
                
        else:
            nwater = nwater +1
            ntot = ntot + 1
        

    elif types[iatom] == 'HW1':

        if types[iatom - 1] == 'OW':

            if z_atom[iatom - 1] - z_ref < 1.5 and z_atom[iatom - 1] - Zmax < 3.0:
                ntot = ntot + 1
                
            else:
                ntot = ntot + 1

        elif types[iatom - 2] == 'OW':

            if z_atom[iatom - 2] - z_ref < 1.5 and z_atom[iatom - 2] - Zmax < 3.0:
                ntot = ntot + 1
                
            else:
                ntot = ntot + 1
                
                
    elif types[iatom] == 'HW2':

        if types[iatom - 1] == 'OW':

            if z_atom[iatom - 1] - z_ref < 1.5 and z_atom[iatom - 1] - Zmax < 3.0:
                ntot = ntot + 1
                
            else:
                ntot = ntot + 1

        elif types[iatom - 2] == 'OW':

            if z_atom[iatom - 2] - z_ref < 1.5 and z_atom[iatom - 2] - Zmax < 3.0:
                ntot = ntot + 1
                
            else:
                ntot = ntot + 1

    elif types[iatom] == 'MW':

        if types[iatom - 3] == 'OW':

            if z_atom[iatom - 3] - z_ref < 1.5 and z_atom[iatom - 3] - Zmax < 3.0:
                ntot = ntot + 1
                
            else:
                ntot = ntot + 1
            

    elif types[iatom][0] == 'C':

        if z_atom[iatom] - z_ref < 1.5 and z_atom[iatom] - Zmax < 3.0:
            ntot = ntot + 1
            nco = nco+1

            if z_ref < z_atom[iatom]:
                z_ref = z_atom[iatom]
                
        else:
            ntot = ntot + 1
            nco = nco+1

    elif types[iatom][0] == 'O':

        if z_atom[iatom - 1] - z_ref < 1.5 and z_atom[iatom - 1] - Zmax < 3.0:
            ntot = ntot + 1
            
        else:
            ntot = ntot + 1

    elif types[iatom][0] == 'D':

        if z_atom[iatom - 2] - z_ref < 1.5 and z_atom[iatom - 2] - Zmax < 3.0:
            ntot = ntot + 1
            
        else:
            ntot = ntot + 1

print(nsilica, nwater, ntot, nco)

line_box = gro_input_file.readline()
box_cut = line_box.split()
x_box = float(box_cut[0])
y_box = float(box_cut[1])
z_box = float(box_cut[2])

gro_input_file.close()

gro_output_file = open(f'../{output}.gro', 'w')

line_title = 'SiO2 usual surface with CO \n'
gro_output_file.write(line_title)

line_blank = '' + str(ntot+3) + '    \n'

gro_output_file.write(line_blank)

count_water = 301
count_carb = 801

for iatom in range(natoms):

    if iatom < 9:
        
        if types[iatom] == 'SZ':
            
            line = '    ' + str(iatom+1) + 'silic' + '   ' + types[iatom] + '     ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) +'\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'SZ1':
            
            line = '    ' + str(iatom+1) + 'silic' + '   ' + types[iatom] + '     ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) +'\n'
            gro_output_file.write(line)
            
    elif iatom < 99:
        
        if types[iatom] == 'SZ':
            
            line = '   ' + str(iatom+1) + 'silic' + '   ' + types[iatom] + '    ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) +'\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'SZ1':
            
            line = '   ' + str(iatom+1) + 'silic' + '   ' + types[iatom] + '    ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) +'\n'
            gro_output_file.write(line)
            

    elif iatom < 999:
        
        if types[iatom] == 'SZ':
            
            line = '  ' + str(iatom+1) + 'silic' + '    ' + types[iatom] + '   ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) +'\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'SZ1':
            
            line = '  ' + str(iatom+1) + 'silic' + '   ' + types[iatom] + '   ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) +'\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'OW':
            
            line = '  ' + str(count_water) + 'water' + '    ' + types[iatom] + '   ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'HW1':
            
            line = '  ' + str(count_water) + 'water' + '   ' + types[iatom] + '   ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'HW2':
            
            line = '  ' + str(count_water) + 'water' + '   ' + types[iatom] + '   ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'MW':
            
            line = '  ' + str(count_water) + 'water' + '    ' + types[iatom] + '   ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            count_water = count_water + 1
            
        elif types[iatom] == 'C':
            
            line = '     ' + str(count_carb) + 'carbo' + '      ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'D':
            
            line = '     ' + str(count_carb) + 'carbo' + '      ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'O':
            
            line = '     ' + str(count_carb) + 'carbo' + '      ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            count_carb = count_carb + 1
            


    elif iatom < 9999:
        
        if types[iatom] == 'OW':
            
            line = ' ' + str(count_water) + 'water' + '    ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'HW1':
            
            line = ' ' + str(count_water) + 'water' + '   ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'HW2':
            
            line = ' ' + str(count_water) + 'water' + '   ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'MW':
            
            line = ' ' + str(count_water) + 'water' + '    ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            count_water = count_water + 1
            
        elif types[iatom] == 'C':
            
            line = ' ' + str(count_carb) + 'carbo' + '     ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
        elif types[iatom] == 'D':
            
            line = ' ' + str(count_carb) + 'carbo' + '     ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            count_carb = count_carb + 1
            
        elif types[iatom] == 'O':
            
            line = ' ' + str(count_carb) + 'carbo' + '     ' + types[iatom] + '  ' + str(
                iatom + 1) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                v_x[iatom]) + "{:8.4f}".format(
                v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
            gro_output_file.write(line)
            
            


    elif iatom < 99999:

        if types[iatom] == 'OW':

            if z_atom[iatom]-z_ref < 1.5 and z_atom[iatom] - Zmax < 3.0:

                line = '' + str(iatom) + 'OW' + '    ' + 'OW' + '' + str(
                    iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                    y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                    v_x[iatom]) + "{:8.4f}".format(
                    v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                gro_output_file.write(line)

        elif types[iatom] == 'HW':

            if types[iatom-1] == 'OW':

                if z_atom[iatom-1] - z_ref < 1.5 and z_atom[iatom-1] - Zmax < 3.0:

                    line = '' + str(iatom) + 'HW1' + '    ' + 'HW1' + '' + str(
                        iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                        y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                        v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                    gro_output_file.write(line)

            elif types[iatom-2] == 'OW':

                if z_atom[iatom-2] - z_ref < 1.5 and z_atom[iatom-2] - Zmax < 3.0:

                    line = '' + str(iatom) + 'HW2' + '    ' + 'HW2' + '' + str(
                        iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                        y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                        v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                    gro_output_file.write(line)

        elif types[iatom] == 'MW':

            if types[iatom-3] == 'OW':

                if z_atom[iatom-3] - z_ref < 1.5 and z_atom[iatom-3] - Zmax < 3.0:

                    line = '' + str(iatom) + 'HW1' + '    ' + 'HW1' + '' + str(
                        iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                        y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(v_x[iatom]) + "{:8.4f}".format(
                        v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                    gro_output_file.write(line)

        elif types[iatom][0] == 'C':

            if z_atom[iatom]-z_ref < 1.5 and z_atom[iatom] - Zmax < 3.0:

                line = '' + str(iatom) + ' C' + '    ' + ' C' + '' + str(
                    iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                    y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                    v_x[iatom]) + "{:8.4f}".format(
                    v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                gro_output_file.write(line)

        elif types[iatom][0] == 'O':

            if z_atom[iatom-1]-z_ref < 1.5 and z_atom[iatom-1] - Zmax < 3.0:

                line = '' + str(iatom) + ' O' + '    ' + ' O' + '' + str(
                    iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                    y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                    v_x[iatom]) + "{:8.4f}".format(
                    v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                gro_output_file.write(line)

        elif types[iatom][0] == 'D':

            if z_atom[iatom-2]-z_ref < 1.5 and z_atom[iatom-2] - Zmax < 3.0:

                line = '' + str(iatom) + ' D' + '    ' + ' D' + '' + str(
                    iatom) + "{:8.3f}".format(x_atom[iatom]) + "{:8.3f}".format(
                    y_atom[iatom]) + "{:8.3f}".format(z_atom[iatom]) + "{:8.4f}".format(
                    v_x[iatom]) + "{:8.4f}".format(
                    v_y[iatom]) + "{:8.4f}".format(v_z[iatom]) + '\n'
                gro_output_file.write(line)

z_altitude = 0.0
initial_temperature = 300.0

initial_x_C = random.random()*x_box + x_ref
initial_y_C = random.random()*y_box + y_ref

for jatom in range(natoms):

    if np.sqrt((x_atom[jatom]-initial_x_C)**2+(y_atom[jatom]-initial_y_C)**2) < 1.5:

        if z_altitude < z_atom[jatom] - Zmax:

            z_altitude = z_atom[jatom] - Zmax

initial_altitude_C = Zmax + z_altitude + 1.5

direction_x_ini_1 = random.uniform(-1,1)
direction_y_ini_1 = random.uniform(-1,1)
direction_z_ini_1 = random.uniform(-1,1)

initial_x_O = initial_x_C + 0.1128*direction_x_ini_1/(np.sqrt(direction_x_ini_1**2+direction_y_ini_1**2+direction_z_ini_1**2))
initial_y_O = initial_y_C + 0.1128*direction_y_ini_1/(np.sqrt(direction_x_ini_1**2+direction_y_ini_1**2+direction_z_ini_1**2))
initial_z_O = initial_altitude_C + 0.1128*direction_z_ini_1/(np.sqrt(direction_x_ini_1**2+direction_y_ini_1**2+direction_z_ini_1**2))

initial_x_D = (12.0*initial_x_C+16.0*initial_x_O)/28.0
initial_y_D = (12.0*initial_y_C+16.0*initial_y_O)/28.0
initial_z_D = (12.0*initial_altitude_C+16.0*initial_z_O)/28.0

locat_vel = 3.72e-2*np.sqrt(initial_temperature)
scale_vel = 2.15e-2*np.sqrt(initial_temperature)

initial_speed = maxwell.rvs(loc = locat_vel, scale = scale_vel, size = 1)

initial_v_x_unscaled = random.uniform(-1,1)
initial_v_y_unscaled = random.uniform(-1,1)
initial_v_z_unscaled = random.uniform(-1,0)

initial_v_x = initial_v_x_unscaled*initial_speed/(np.sqrt(initial_v_x_unscaled**2+initial_v_y_unscaled**2+initial_v_z_unscaled**2))
initial_v_y = initial_v_y_unscaled*initial_speed/(np.sqrt(initial_v_x_unscaled**2+initial_v_y_unscaled**2+initial_v_z_unscaled**2))
initial_v_z = initial_v_z_unscaled*initial_speed/(np.sqrt(initial_v_x_unscaled**2+initial_v_y_unscaled**2+initial_v_z_unscaled**2))

#initial_v_x = initial_speed*0.0
#initial_v_y = initial_speed*0.0
#initial_v_z = -initial_speed


line = ' '+str(count_carb) + 'carbo' + '    ' + ' C' + '  ' + str(
            natoms + 1) + "{:8.3f}".format(initial_x_C) + "{:8.3f}".format(
    initial_y_C) + "{:8.3f}".format(initial_altitude_C)+ "{:8.4f}".format(initial_v_x[0]) + "{:8.4f}".format(
    initial_v_y[0]) + "{:8.4f}".format(initial_v_z[0]) +'\n'
gro_output_file.write(line)

line = ' '+str(count_carb) + 'carbo' + '   ' + '  O' + '  ' + str(
            natoms + 2) + "{:8.3f}".format(initial_x_O) + "{:8.3f}".format(
    initial_y_O) + "{:8.3f}".format(initial_z_O)+ "{:8.4f}".format(initial_v_x[0]) + "{:8.4f}".format(
    initial_v_y[0]) + "{:8.4f}".format(initial_v_z[0]) +'\n'
gro_output_file.write(line)

line = ' '+str(count_carb) + 'carbo' + '   ' + '  D' + '  ' + str(
            natoms + 3) + "{:8.3f}".format(initial_x_D) +"{:8.3f}".format(
    initial_y_D) + "{:8.3f}".format(initial_z_D)+ "{:8.4f}".format(initial_v_x[0]) + "{:8.4f}".format(
    initial_v_y[0]) + "{:8.4f}".format(initial_v_z[0]) +'\n'
gro_output_file.write(line)

line = '   '+"{0:.5f}".format(x_box)+'   '+"{0:.5f}".format(y_box)+'   '+"{0:.5f}".format(z_box)+'   '
gro_output_file.write(line)

gro_output_file.close()


