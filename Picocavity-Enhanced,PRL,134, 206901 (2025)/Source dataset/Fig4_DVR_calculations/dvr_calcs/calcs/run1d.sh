#!/bin/bash


echo "---- HYDROGEN ----"
echo "HIGH density"
dens="high"
for dtip in 11.4 9.4 7.9 7.4 ; do
    onedim_model.py --isotope=H --density=${dens} ${dtip} dvr1d_h.json
done

echo ""
echo "LOW density"
dens="low"
for dtip in 11.4 9.4 7.9 7.4 ; do
    onedim_model.py --isotope=H --density=${dens} ${dtip} dvr1d_h.json
done


echo "---- DEUTERIUM ----"
echo "HIGH density"
dens="high"
for dtip in 11.4 9.4 7.9 7.4 ; do
    onedim_model.py --isotope=D --density=${dens} ${dtip} dvr1d_d.json
done

echo ""
echo "LOW density"
dens="low"
for dtip in 11.4 9.4 7.9 7.4 ; do
    onedim_model.py --isotope=D --density=${dens} ${dtip} dvr1d_d.json
done
