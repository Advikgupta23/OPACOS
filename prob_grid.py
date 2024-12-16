import numpy as np
import pandas as pd
#from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages 
import matplotlib.pyplot as plt
from isobuild import func
from multiprocessing import Pool
import Extinspiral as exbv_spiral
import run_extinction as exbv

#from constants import age, feh, nstars, imf_slope
# define array of [Fe/H] from -3 to 0.5 in steps of 0.1
feh = np.arange(-3, 0.6, 0.1)
#feh = np.array([-3.])  # smaller range of [Fe/H]
# define array of ages from 0.5 to 13 in steps of 0.5. Units of Gyr
age = np.arange(0.5, 14.5, 0.125)
#age = np.array([11.1])  # smaller range of ages

nstars = 100000  # number of stars
imf_slope = -2.35  # IMF slope. -2.35 is Salpeter pure power-law
imf_type = 'kroupa' # You can use three IMF systems: * salpeter 
                      #                                * kroupa  
                      #                                * chabrierlognormal

Phot_sys = '2mass' # You can use * 'GAIA_EDR3' or * '2mass'
                             
fai_lim = 14.0  # faint magnitude limit 20.0 19.0
bri_lim = 12.0   # bright magnitude limit 3.0 12.0 
col_slo = 0.0  # colour slope 0.20
blu_lim = -100.0  # lower colour limit -2.5 0
red_lim = 100.0  # upper colour limit 5.1 1.5

low_l = -9  # log(g) cut. E.g. between 1 and 4. If you write
hig_l = -9  # -9. and -9, then log(g) cut is ignored.

midi = 10.0  # minimum distance
madi = 10000.0  # maximum distance
dedi = 100.0  # distance step

extinction_mode = 'in_plane' # 'in_pkane' for the computation of extinction using extinction due to spiral arms and schlegel maps
                             # 'out_plane' for the computation of extinction using extinction using general schlegel maps
R_BP = 3.1
R_RP = 2.2
R_G = 2.7

R_J = 0.72
R_K = 0.306
R_V = 3.1

#G_long = 309.1
#G_lat = 14.97

G_long = 348.0
G_lat = 18.0

def read_columns(filepath):
    
    chunk = pd.read_csv(filepath, sep = '\s+', comment='#', header=None,
                         names = ['mini','mfin','age','feh','distance','l','b','px','py','pz','popid','exbv_schlegel'],
                         usecols = [0,1,2,3,4,5,6,7,8,9,10,11], chunksize = 200)
    data = pd.concat(chunk)
    
    return data

if extinction_mode == 'in_plane':
    ex_bv = exbv_spiral.list_ebv(G_long,G_lat)
elif extinction_mode == 'out_plane':
    ex_bv = exbv.list_ebv(G_long,G_lat)
# print(ex_bv)
# exit()

#filename = 'galaxia.dat'
#data = read_columns(filename)

#l = data['l']
#b = data['b']
#data.sort_values(by=['distance'], ascending=[True], inplace=True, ignore_index=True)

#ex_bv = data['exbv_schlegel']
#distance_model = data['distance']

def inter_rout(x1, y1, x2, y2, x): 
    y = ( (y2-y1)/(x2-x1) ) * (x-x1) + y1
    return y

with open('test_results.dat', 'w') as f:

    f.write(f'# stars {nstars}\n')
    if imf_type == 'salpeter':
       f.write(f'# Salpeter IMF = {imf_slope:.4f}\n')
    else:
       f.write(f'# IMF = {imf_type}\n')
    f.write(f'# Photometric System used:{Phot_sys}\n')
    f.write(f'# faint magnitude limit = {fai_lim:.4f}\n')
    f.write(f'# bright magnitude limit = {bri_lim:.4f}\n')
    f.write(f'# colour slope = {col_slo:.4f}\n')
    f.write(f'# blue cut = {blu_lim:.4f} and red cut = {red_lim:.4f}\n')
    f.write(f'# sampled distances from {midi:.1f} to {madi:.1f} in steps of {dedi:.1f} pc\n')
    if low_l < -8 and hig_l < -8:
       f.write('# no logg cut\n')
    else:
       f.write(f'# {low_l:.2f} < logg < {hig_l:.2f}\n')
    f.write('#------------------------------------------------------\n')
    f.write('#                                                      \n')
    f.write('#    Age(Gyr)       [Fe/H]    distance     Prob        \n')

    
    #def save_multi_image(filename):
    #    pp = PdfPages(filename)
    #    fig_nums = plt.get_fignums()
    #    figs = [plt.figure(n) for n in fig_nums]
    #    for fig in figs:
    #        fig.savefig(pp, format='pdf')
    #    pp.close()    
    
    
    Func = func(age, feh, nstars, imf_slope, imf_type)

    # Perform any necessary setup or initialization
    Func.prepare_computations()
    
    # Computation cycle through age-metallicity grid
    #for distance in np.arange(midi, madi, dedi):
        #for index in range(len(distance_model)):
            #          if distance_model[index] > dis:
            #             i = index
            #             break	
            #ex_bv_dist = inter_rout(distance_model[i-1],ex_bv[i-1],distance_model[i],ex_bv[i],dis)
            #ex_bv_long = inter_rout(l[i-1],ex_bv[i-1],l[i],ex_bv[i],G_long)
            #ex_bv_lat = inter_rout(b[i-1],ex_bv[i-1],b[i],ex_bv[i],G_lat)
            #extinc_bv = np.sum(ex_bv_dist + ex_bv_long + ex_bv_lat)/3
            #extinc_bv = round(extinc_bv,4)
        
            #extinc_bv = np.sum(ex_bv[i-1] + ex_bv[i])/2
            #extinc_bv = round(extinc_bv,4)
        
            #blu_limit = blu_lim + extinc_bv*(R_BP-R_RP)
            #blu_limit = round(blu_limit,4)
            #red_limit = red_lim + extinc_bv*(R_BP-R_RP)
            #red_limit = round(red_limit,4)
            #fai_limit = fai_lim + extinc_bv*(R_G)
            #fai_limit = round(fai_limit,4)
            #bri_limit = bri_lim + extinc_bv*(R_G)
            #bri_limit = round(bri_limit,4) 
    
    for tau in age:
        
        for metal in feh:
           
           #pool = Pool(processes=10)
           #lst = [(1,2),(3,4),(5,6)]
           color, mag, logg, nstars = Func.compute_parameters(tau, metal, Phot_sys)
           
           #plt.figure()
           #plt.scatter(color,mag,s=0.5)
           #plt.xlabel('J-K')
           #plt.ylabel(r'$M_{V}$')
           #plt.gca().invert_yaxis()
           
           
           #plt.show()
           #fig = plt.gcf()
           #plt.savefig(f'Plots/isochrone_age({tau})_metal({metal}).png')
           count = 0
           for distance in np.arange(midi, madi, dedi):
               
               dis = float(distance)
               apparent = mag + 5.0 * np.log10(dis) - 5.0

               # plt.figure()
               # plt.scatter(color,apparent,s=0.5)
               # plt.xlabel('J-K')
               # plt.ylabel(r'$m_{V}$')
               # plt.gca().invert_yaxis()
               # plt.style.use('dark_background')

               # #plt.show()
               # fig = plt.gcf()
               # plt.savefig(f'Plots/isochrone_age({tau})_metal({metal})_distance({dis}).png')
               # plt.close()

               if Phot_sys == '2mass': 
                  
                  extinc_bv = ex_bv[count]
                  count+=1

                  blu_limit = blu_lim + extinc_bv*(R_J-R_K)
                  blu_limit = round(blu_limit,4)
                  red_limit = red_lim + extinc_bv*(R_J-R_K)
                  red_limit = round(red_limit,4)
                  fai_limit = fai_lim + extinc_bv*(R_V)
                  fai_limit = round(fai_limit,4)
                  bri_limit = bri_lim + extinc_bv*(R_V)
                  bri_limit = round(bri_limit,4)
                  
                  if low_l < -8 and hig_l < -8:
                     passed = np.sum((color >= blu_limit) & (color <= red_limit) &
                                    (apparent <= fai_limit - col_slo * color) &
                                    (apparent >= bri_limit - col_slo * color))
                  else:
                     passed = np.sum((color >= blu_limit) & (color <= red_limit) &
                                    (apparent <= fai_limit - col_slo * color) &
                                    (apparent >= bri_limit - col_slo * color) &
                                    (logg >= low_l) & (logg <= hig_l))
                
                  prob = float(passed) / float(nstars)
           
                  f.write(f'       {tau:.3f}        {metal:.2f}      {distance:.2f}        {prob:.6f}\n')

               if Phot_sys == 'GAIA_EDR3':
                  #for index in range(len(distance_model)):
                  #    if distance_model[index] > dis:
                  #       i = index
                  #       break	
                  #ex_bv_dist = inter_rout(distance_model[i-1],ex_bv[i-1],distance_model[i],ex_bv[i],dis)
                  #ex_bv_long = inter_rout(l[i-1],ex_bv[i-1],l[i],ex_bv[i],G_long)
                  #ex_bv_lat = inter_rout(b[i-1],ex_bv[i-1],b[i],ex_bv[i],G_lat)
                  #extinc_bv = np.sum(ex_bv_dist + ex_bv_long + ex_bv_lat)/3
                  #extinc_bv = round(extinc_bv,4)
                  
                  #extinc_bv = np.sum(ex_bv[i-1] + ex_bv[i])/2
                  
                  extinc_bv = ex_bv[count]
                  count+=1
                  
                  blu_limit = blu_lim + extinc_bv*(R_BP-R_RP)
                  blu_limit = round(blu_limit,4)
                  red_limit = red_lim + extinc_bv*(R_BP-R_RP)
                  red_limit = round(red_limit,4)
                  fai_limit = fai_lim + extinc_bv*(R_G)
                  fai_limit = round(fai_limit,4)
                  bri_limit = bri_lim + extinc_bv*(R_G)
                  bri_limit = round(bri_limit,4)
                  #print(extinc_bv)
                                                     
                  if low_l < -8 and hig_l < -8:
                     
                     passed = np.sum((color >= blu_limit) & (color <= red_limit) &
                                    (apparent <= fai_limit - col_slo * color) &
                                    (apparent >= bri_limit - col_slo * color))
                  else:
                     passed = np.sum((color >= blu_limit) & (color <= red_limit) &
                                    (apparent <= fai_limit - col_slo * color) &
                                    (apparent >= bri_limit - col_slo * color) &
                                    (logg >= low_l) & (logg <= hig_l))
                
                  prob = float(passed) / float(nstars)
           
                  f.write(f'       {tau:.3f}        {metal:.2f}      {distance:.2f}        {prob:.6f}\n')                 
          
           
           
        print(f'done with grid at {tau} Gyr')
    #filename = "plots.pdf"
    #save_multi_image(filename)cd 
