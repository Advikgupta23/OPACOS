#!/bin/bash

#x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x

##################################################################
##################################################################
######################## INPUT ARGUMENTS #########################
##################################################################
##################################################################

galactic_longitude=67.0
galactic_latitude=-13.0
survey_area=10.0
Phot_sys='2mass'  # Right now you can use two type of photometric system (passed in data_cube.py) : * 'GAIA_EDR3'
                  #                                                                                 * '2mass'.
                  # Our stars of interest are part of GALAH survey with selection function in 2mass filters.

imf_type='kroupa' # You can use three IMF systems while generating the probability (passed in data_cube.py) : * 'salpeter' 
                    #                                                                                           * 'kroupa'  
                    #                                                                                            * 'chabrierlognormal'.
extinction_mode='in_plane' # 'in_plane' for the computation of extinction using extinction due to spiral arms and schlegel maps (passed in data_cube.py).
                            # 'out_plane' for the computation of extinction using extinction using general schlegel maps (passed in data_cube.py).
age_distribution=False #If you want to infer model age distribution for the section of targeted sky set it to be True otherwise False. 
feh_distribution=True  #If you want to infer model feh distribution for the section of targeted sky set it to be True otherwise False.

#x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x


python3 modify_myparameterfile.py "$galactic_longitude" "$galactic_latitude" "$survey_area" && echo -e "\n\n***Successfully modified the parameter file for GALAXIA simulation***\n\n"

galaxia -r ../GalaxiaData/Examples/myparameterfile && python3 galaxy1.py
 
echo -e "\n\n***Successful GALAXIA simulation; model saved to galaxia.dat***\n\n"

time python3 data_cube.py "$galactic_longitude" "$galactic_latitude" "$Phot_sys" "$imf_type" "$extinction_mode" && echo -e "\n\n***Successfully compiled the probability grid file***\n\n"

time python3 PDF_estimation.py "$galactic_longitude" "$galactic_latitude" "$age_distribution" "$feh_distribution" && echo -e "\n\n***The model, data distributions, and associated data have been successfully saved in the RESULTS folder.***\n\n"
