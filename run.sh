#!/bin/bash

galactic_longitude=69.0
galactic_latitude=-38.0

python3 modify_myparameterfile.py "$galactic_longitude" "$galactic_latitude" && echo -e "\n\n***Successfully modified the parameter file for GALAXIA simulation***\n\n"

galaxia -r ../GalaxiaData/Examples/myparameterfile && python3 galaxy1.py
 
echo -e "\n\n***Successful GALAXIA simulation; model saved to galaxia.dat***\n\n"

#python3 data_cube.py "$galactic_longitude" "$galactic_latitude" && echo -e "\n\n***Successfully compiled the probability grid file***\n\n"

python3 PDF_estimation.py "$galactic_longitude" "$galactic_latitude" && echo -e "\n\n***The model, data distributions, and associated data have been successfully saved in the RESULTS folder.***\n\n"
