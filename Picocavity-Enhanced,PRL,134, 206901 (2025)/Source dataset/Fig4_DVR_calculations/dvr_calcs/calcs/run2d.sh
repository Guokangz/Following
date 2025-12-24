#!/bin/bash


echo "---- HYDROGEN ----"
echo "HIGH density"
dens="high"
for dtip in 9.4 7.4 ; do
    twodim_model.py --isotope=H --density=${dens} ${dtip} dvr2d_h.json
done


echo "---- DEUTERIUM ----"
echo "HIGH density"
dens="high"
for dtip in 9.4 7.4 ; do
    twodim_model.py --isotope=D --density=${dens} ${dtip} dvr2d_d.json
done

