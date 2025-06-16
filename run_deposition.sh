#! /usr/bin/bash

python3 adding_one_water_new.py

gmx grompp -maxwarn 3 -f relax_20K_quick.mdp -c silic_surface_3H2O.gro -p topol_3H2O.top -o quick_silic_surface_3H2O_relax_20K.tpr

gmx mdrun -s quick_silic_surface_3H2O_relax_20K.tpr -deffnm quick_silic_surface_3H2O_relax_20K
