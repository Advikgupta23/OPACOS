import numpy as np
import ebf

data = ebf.read('/home/advik/GalaxiaData/Examples/galaxy1.ebf','/')
rad = data['rad']
age = data['age']
age = pow(10,data['age'])/pow(10,9)
feh = data['feh']
glon = data['glon']
glat = data['glat']
mass_fin = data['mact']
mass_in = data['smass']
px = data['px']
py = data['py']
pz = data['pz']
popid = data['popid']
exbv_schlegel = data['exbv_schlegel']
log_g = data['grav']

with open('galaxia.dat', 'w') as f:
     
     f.write(f'#Mass_initial       Mass_final    Age(Gyr)       [Fe/H]      distance(pc)     Galactic_long.     Galactic_lat.        px               py              pz            popid          exbv_schlegel        log_g \n')
     for i in range(len(data['rad'])):
         
         f.write(f'    {mass_in[i]:.4f}            {mass_fin[i]:.4f}       {age[i]:.4f}        {feh[i]:.4f}       {1000*rad[i]:.4f}           {glon[i]:.4f}           {glat[i]:.4f}         {px[i]:.4f}          {py[i]:.4f}          {pz[i]:.4f}          {popid[i]:.4f}             {exbv_schlegel[i]:.4f}             {log_g[i]:.4f} \n')
         
         
         
         
