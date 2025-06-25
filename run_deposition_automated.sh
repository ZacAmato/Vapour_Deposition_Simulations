#! /usr/bin/bash

COUNTER=1

python3 adding_one_water_new.py

printf $COUNTER

gmx grompp -maxwarn 3 -f relax_20K_quick.mdp -c quick_silic_surface_${COUNTER}H2O_relax_20K.gro -p topol_H2O.top -o quick_silic_surface_${COUNTER}H2O_relax_20K.tpr

gmx mdrun -s quick_silic_surface_${COUNTER}H2O_relax_20K.tpr -deffnm quick_silic_surface_${COUNTER}H2O_relax_20K

# Add python script to check water attached

# Then run again for another 5 ps or less

# if its not attached then do another 



#########

sed -i "7 s/quick_silic_surface_relax_20K/quick_silic_surface_${COUNTER}H2O_relax_20K/" "adding_one_water_new.py"

let Nwater=$COUNTER+1

printf $Nwater

sed -i "9 s/quick_silic_surface_${COUNTER}H2O_relax_20K/quick_silic_surface_${Nwater}H2O_relax_20K/" "adding_one_water_new.py"

python3 adding_one_water_new.py

sed -i "77 s/water ${COUNTER}/water ${Nwater}/" "topol_H2O.top"

gmx grompp -maxwarn 3 -f relax_20K_quick.mdp -c quick_silic_surface_${Nwater}H2O_relax_20K.gro -p topol_H2O.top -o quick_silic_surface_${Nwater}H2O_relax_20K.tpr

gmx mdrun -s quick_silic_surface_${Nwater}H2O_relax_20K.tpr -deffnm quick_silic_surface_${Nwater}H2O_relax_20K

##########

