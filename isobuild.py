import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages 
import matplotlib.pyplot as plt
import os   
import multiprocessing as mp
import time
import imf
#from test import age, feh, nstars, imf_slope

#def compute_inter_rout(mge1, val1, mge2, val2, tau):
#    return inter_rout(mge1, val1, mge2, val2, tau)
                      
#def compute_inter_rout_b(bmet1, val1, bmet2, val2, metal):
#    return inter_rout(bmet1, val1, bmet2, val2, metal) 

#def inter_rout(x1, y1, x2, y2, x):
#    f = interp1d([x1, x2], [y1, y2], fill_value="extrapolate")
#    return f(x)

class func:

        def __init__(self, age, feh, nstars, imf_slope, imf_type):
	    	    
            self.age = age
            self.feh = feh
            self.nstars = nstars
            self.imf_slope = imf_slope
            self.imf_type = imf_type
            
        def prepare_computations(self):
            
            pass
            
        def compute_parameters(self, tau, metal, Phot_sys):
            
           
            color, mag, logg, nstars = self.Isobuild(tau, metal, Phot_sys)
            return color, mag, logg, nstars
           
        def Isobuild(self, tau, metal, Phot_sys):

            
            def read_isochrone_models(filename1):
                bmod = []  # List to store the names of isochrone models
                bmet = []  # List to store the corresponding metallicity values
	    
                with open(filename1, 'r') as file:
                     for line in file:
                         line = line.strip()  # Remove leading/trailing whitespace
                         if line:  # Skip empty lines
                            parts = line.split()  # Split the line into parts
                            if len(parts) == 2:  # Ensure there are exactly two parts
                               bmod.append(parts[0])  # Add the model name to the list
                               bmet.append(float(parts[1]))  # Add the metallicity value to the list
                return bmod, bmet
	    
            
            if Phot_sys == '2mass':
               filename1 = 'set.input.basti'
	    
            if Phot_sys == 'GAIA_EDR3':
               filename1 = 'set.input.basti.GAIA_EDR3'
            
            bmod, bmet = read_isochrone_models(filename1)
		
		
	    # path to where isochone models are located
            path   = 'ncan_eta00/'
            logsun = 4.44           # Solar log(g)
            Tsun   = 5777.          # Solar Teff
	
            order  = 'agemet'   # agemet or metage for isochrone interpolation order
	     		    # the two method should be equivalent
	
			   	
	    # minimum tau in BASTI is usually 30 Myr (0.03 Gry), but for the most metal poor
	    # models minimum tau might be 40 Myr or 90 Myr. Maximum tau in BASTI is 15 Gyr
	        
            if tau < 0.03 :
               tau = 0.0301
            if tau < 0.04 and metal < -2.489:
               tau = 0.0401
            if tau < 0.09 and metal < -2.789:
               tau = 0.0901
            if tau > 15.:
               tau = 14.998
	    
	    # BASTI isochrones are limited to -3.27 < [Fe/H] < +0.5.
            if Phot_sys == '2mass':
    
               if metal < -3.27:
                  metal = -3.269
    
               if metal > 0.5:
                  metal = 0.5
    
            if Phot_sys == 'GAIA_EDR3':
       
               if metal < -3.2:
                  metal = -3.19
       
               if metal > 0.45:
                  metal = 0.45
		
		
            j = 0
            for ff in range(len(bmod)):  # Equivalent to IDL's for loop range 0,22
                if metal <= bmet[ff] and metal > bmet[ff + 1]:
                   j = ff	
	
            def sort_filenames(filenames):
            # Extract numeric part and create a list of tuples (numeric_part, original_filename)
                filenames_with_numeric = [(int(f.split('z')[0]), f) for f in filenames]

                # Sort the list of tuples by the numeric part
                sorted_filenames_with_numeric = sorted(filenames_with_numeric)
    
                # Extract the sorted filenames
                sorted_filenames = [f[1] for f in sorted_filenames_with_numeric]
                
                return sorted_filenames
               		
	    
	    # Bracket in tau at bmet[j]
            if Phot_sys == '2mass':
            	
               directory1 = os.path.join(path, bmod[j], '2mass')
               listA = [filename for filename in os.listdir(directory1) if filename.startswith('w') and filename.endswith('_c03hbs.gz')]
               listA.sort()
            
               for filename in listA:
                   
                   filepath = os.path.join(directory1, filename)
                   junk, eta = filepath[:41], filepath[41:46]
                   eta = float(eta.strip()) / 1000.0
                   if eta > tau:
                      break  
            
            if Phot_sys == 'GAIA_EDR3':
               
               directory1 = os.path.join(path, bmod[j])
               listA = [filename for filename in os.listdir(directory1) if  filename.endswith('D0E0.isc_gaia-dr3')]
               listA = sort_filenames(listA)

               for filename in listA:
                   filepath = os.path.join(directory1, filename)
                   junk, eta_temp = filepath[:19], filepath[19:25]
                   eta = int(eta_temp.partition("z")[0]) / 1000.0

                   if eta > tau:
                      break
            #print(filepath)
            #print(eta, eta_temp)
	    #ncan_eta00/fehm327/2mass/wz105y245soe0.t600090_c03hbs.gz
	    #ncan_eta00/FEHm320/30z0000100y247P00O1D0E0.isc_gaia-dr3
	    
	    # ExitA
            modA2 = filename
            mgeA2 = eta
            
	    # Since listA is a list of filenames, we can find the index of modA2 in listA
            index_modA2 = listA.index(modA2)
            modA1 = listA[index_modA2 - 1]
          	
            filepath = os.path.join(directory1, modA1)
            
            if Phot_sys == '2mass':
             
               junk, eta = filepath[:41], filepath[41:46]
               eta = float(eta.strip()) / 1000.0
	    
            if Phot_sys == 'GAIA_EDR3':
	       
	       
               junk, eta_temp = filepath[:19], filepath[19:25]
               eta = int(eta_temp.partition("z")[0]) / 1000.0
	       
            mgeA1 = eta
	
            # braket in tau at bmet[j+1]
            if Phot_sys == '2mass':
            	
               directory2 = os.path.join(path, bmod[j+1], '2mass')
               listB = [filename2 for filename2 in os.listdir(directory2) if filename2.startswith('w') and filename2.endswith('_c03hbs.gz')]
               listB.sort()
            
            if Phot_sys == 'GAIA_EDR3':
               
               directory2 = os.path.join(path, bmod[j+1])
               listB = [filename2 for filename2 in os.listdir(directory2) if  filename2.endswith('D0E0.isc_gaia-dr3')]
               listB = sort_filenames(listB)	
               
            for filename2 in listB:
                filepath = os.path.join(directory2, filename2)
		    
                if Phot_sys == '2mass':
                   junk, eta = filepath[:41], filepath[41:46]
                   eta = float(eta.strip()) / 1000.0
                   if eta > tau:
                      break
                
                if Phot_sys == 'GAIA_EDR3':
                   
                   junk, eta_temp = filepath[:19], filepath[19:25]
                   eta = int(eta_temp.partition("z")[0]) / 1000.0                
                   if eta > tau:
                      break
	
	    # ExitB
            modB2 = filename2
            mgeB2 = eta
	
	    # Since listB is a list of filenames, we can find the index of modB2 in listB
            index_modB2 = listB.index(modB2)
            modB1 = listB[index_modB2 - 1]
	
	
            filepath = os.path.join(directory2, modB1)
            
            if Phot_sys == '2mass':
             
               junk, eta = filepath[:41], filepath[41:46]
               eta = float(eta.strip()) / 1000.0
	    
            if Phot_sys == 'GAIA_EDR3':
	       
               junk, eta_temp = filepath[:19], filepath[19:25]
               eta = int(eta_temp.partition("z")[0]) / 1000.0	        
          
       
            mgeB1 = eta
        
            def inter_rout(x1, y1, x2, y2, x): 
                y = ( (y2-y1)/(x2-x1) ) * (x-x1) + y1
                return y

            
	    # Define the function for reading the isochrone files
            if Phot_sys == '2mass':
               def read_iso_file(file_path):
	           # Read the file into a pandas DataFrame
                   chunk = pd.read_csv(file_path, sep='\s+', comment='#', header=None,
                                      names=['min', 'mfin', 'lum', 'logT', 'Mv', 'ub', 'bv', 'vi', 'vr', 'vj', 'vk'],
                                      usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],chunksize=200)
                   data = pd.concat(chunk)
                   return data
               MODA1 = os.path.join(directory1, modA1)
               MODA2 = os.path.join(directory1, modA2)
               MODB1 = os.path.join(directory2, modB1)
               MODB2 = os.path.join(directory2, modB2)
	
               
               dataA1 = read_iso_file(MODA1)
               dataA2 = read_iso_file(MODA2)
               dataB1 = read_iso_file(MODB1)
               dataB2 = read_iso_file(MODB2)

               # Calculate Mk and jk for each dataset
               
              # dataA1['Mk'] = dataA1['Mv'] - dataA1['vk']
              # dataA2['Mk'] = dataA2['Mv'] - dataA2['vk']
              # dataB1['Mk'] = dataB1['Mv'] - dataB1['vk']
              # dataB2['Mk'] = dataB2['Mv'] - dataB2['vk']
               
               dataA1['jk'] = dataA1['vk'] - dataA1['vj']
               dataA2['jk'] = dataA2['vk'] - dataA2['vj']
               dataB1['jk'] = dataB1['vk'] - dataB1['vj']
               dataB2['jk'] = dataB2['vk'] - dataB2['vj']
               
	       # Check for matching errors in the number of elements
               if len(dataA1) != len(dataA2):
                  print('MATCHING ERROR: dataA1 and dataA2 have different lengths')
               if len(dataB1) != len(dataB2):
                  print('MATCHING ERROR: dataB1 and dataB2 have different lengths')
               if len(dataA1) != len(dataB2):
                  print('MATCHING ERROR: dataA1 and dataB2 have different lengths')
	
	       # Define the number of elements in the isochrone
               nbasti = len(dataA1)
	
	       # Initialize arrays for storing interpolated values
               mini = np.zeros(nbasti)
               mfin = np.zeros(nbasti)
               lum = np.zeros(nbasti)
               logT = np.zeros(nbasti)
               Mv = np.zeros(nbasti)
               jk = np.zeros(nbasti)
               
               #def worker(queue, mgeA1, dataA1_min, mgeA2, dataA2_min, tau):
               #    result = inter_rout(mgeA1, dataA1_min, mgeA2, dataA2_min, tau)
               #    queue.put(result)
               # Create a queue for communication between processes
               #queue = Queue()              
               start = time.time()
               for wc in range(nbasti):
                   if order == 'agemet':
                      
                      #i_minA = Process(target=worker, args=(queue, mgeA1, dataA1['min'].iloc[wc], mgeA2, dataA2['min'].iloc[wc], tau))
                      #i_mfinA = Process(target=worker, args=(queue, mgeA1, dataA1['mfin'].iloc[wc], mgeA2, dataA2['mfin'].iloc[wc], tau))
                      #i_lumA = Process(target=worker, args=(queue, mgeA1, dataA1['lum'].iloc[wc], mgeA2, dataA2['lum'].iloc[wc], tau))
                      #i_logTA = Process(target=worker, args=(queue, mgeA1, dataA1['logT'].iloc[wc], mgeA2, dataA2['logT'].iloc[wc], tau))
                      #i_MkA = Process(target=worker, args=(queue, mgeA1, dataA1['Mk'].iloc[wc], mgeA2, dataA2['Mk'].iloc[wc], tau))
                      #i_jkA = Process(target=worker, args=(queue, mgeA1, dataA1['jk'].iloc[wc], mgeA2, dataA2['jk'].iloc[wc], tau))
                      #i_minA.start()
                      #i_mfinA.start()
                      #i_lumA.start()
                      #i_logTA.start()
                      #i_MkA.start()
                      #i_jkA.start()
                      # Wait for the process to finish and get the result
                      #i_minA.join()
                      #i_minA = queue.get()                     
                      #i_mfinA.join()
                      #i_mfinA = queue.get()   
                      #i_lumA.join()
                      #i_lumA = queue.get()   
                      #i_logTA.join()
                      #i_logTA = queue.get()   
                      #i_MkA.join()
                      #i_MkA = queue.get()   
                      #i_jkA.join()
                      #i_jkA = queue.get()   


                      
                      # Create a pool of processes
                      #pool = mp.Pool(mp.cpu_count())
                      
                      # Compute i_minA, i_mfinA, i_lumA, etc. concurrently
                      #i_minA = pool.apply_async(compute_inter_rout, args=(mgeA1, dataA1['min'].iloc[wc], mgeA2, dataA2['min'].iloc[wc], tau))
                      #i_mfinA = pool.apply_async(compute_inter_rout, args=(mgeA1, dataA1['mfin'].iloc[wc], mgeA2, dataA2['mfin'].iloc[wc], tau))
                      #i_lumA = pool.apply_async(compute_inter_rout, args=(mgeA1, dataA1['lum'].iloc[wc], mgeA2, dataA2['lum'].iloc[wc], tau))
                      #i_logTA = pool.apply_async(compute_inter_rout, args=(mgeA1, dataA1['logT'].iloc[wc], mgeA2, dataA2['logT'].iloc[wc], tau))
                      #i_MkA = pool.apply_async(compute_inter_rout, args=(mgeA1, dataA1['Mk'].iloc[wc], mgeA2, dataA2['Mk'].iloc[wc], tau))
                      #i_jkA = pool.apply_async(compute_inter_rout, args=(mgeA1, dataA1['jk'].iloc[wc], mgeA2, dataA2['jk'].iloc[wc], tau))
                      
                      #i_minB = pool.apply_async(compute_inter_rout, args=(mgeB1, dataB1['min'].iloc[wc], mgeB2, dataB2['min'].iloc[wc], tau))
                      #i_mfinB = pool.apply_async(compute_inter_rout, args=(mgeB1, dataB1['mfin'].iloc[wc], mgeB2, dataB2['mfin'].iloc[wc], tau))
                      #i_lumB = pool.apply_async(compute_inter_rout, args=(mgeB1, dataB1['lum'].iloc[wc], mgeB2, dataB2['lum'].iloc[wc], tau))
                      #i_logTB = pool.apply_async(compute_inter_rout, args=(mgeB1, dataB1['logT'].iloc[wc], mgeB2, dataB2['logT'].iloc[wc], tau))
                      #i_MkB = pool.apply_async(compute_inter_rout, args=(mgeB1, dataB1['Mk'].iloc[wc], mgeB2, dataB2['Mk'].iloc[wc], tau))
                      #i_jkB = pool.apply_async(compute_inter_rout, args=(mgeB1, dataB1['jk'].iloc[wc], mgeB2, dataB2['jk'].iloc[wc], tau))
                      
                      # Wait for all computations to finish
                      #i_minA = i_minA.get()
                      #_mfinA = i_mfinA.get()
                      #i_lumA = i_lumA.get()
                      #i_logTA = i_logTA.get()
                      #i_MkA = i_MkA.get()
                      #i_jkA = i_jkA.get()
                      
                      #i_minB = i_minB.get()
                      #i_mfinB = i_mfinB.get()
                      #i_lumB = i_lumB.get()
                      #i_logTB = i_logTB.get()
                      #i_MkB = i_MkB.get()
                      #i_jkB = i_jkB.get()
                      
                      # Perform remaining calculations
                      #i_min = compute_inter_rout_b(bmet[j], i_minA, bmet[j+1], i_minB, metal)
                      #i_mfin = compute_inter_rout_b(bmet[j], i_mfinA, bmet[j+1], i_mfinB, metal)
                      #i_lum = compute_inter_rout_b(bmet[j], i_lumA, bmet[j+1], i_lumB, metal)
                      #i_logT = compute_inter_rout_b(bmet[j], i_logTA, bmet[j+1], i_logTB, metal)
                      #i_Mk = compute_inter_rout_b(bmet[j], i_MkA, bmet[j+1], i_MkB, metal)
                      #i_jk = compute_inter_rout_b(bmet[j], i_jkA, bmet[j+1], i_jkB, metal)
                      
                      i_minA = inter_rout(mgeA1, dataA1['min'].iloc[wc], mgeA2, dataA2['min'].iloc[wc], tau)
                      i_mfinA = inter_rout(mgeA1, dataA1['mfin'].iloc[wc], mgeA2, dataA2['mfin'].iloc[wc], tau)
                      i_lumA = inter_rout(mgeA1, dataA1['lum'].iloc[wc], mgeA2, dataA2['lum'].iloc[wc], tau)
                      i_logTA = inter_rout(mgeA1, dataA1['logT'].iloc[wc], mgeA2, dataA2['logT'].iloc[wc], tau)
                      i_MvA = inter_rout(mgeA1, dataA1['Mv'].iloc[wc], mgeA2, dataA2['Mv'].iloc[wc], tau)
                      i_jkA = inter_rout(mgeA1, dataA1['jk'].iloc[wc], mgeA2, dataA2['jk'].iloc[wc], tau)

                      i_minB = inter_rout(mgeB1, dataB1['min'].iloc[wc], mgeB2, dataB2['min'].iloc[wc], tau)
                      i_mfinB = inter_rout(mgeB1, dataB1['mfin'].iloc[wc], mgeB2, dataB2['mfin'].iloc[wc], tau)
                      i_lumB = inter_rout(mgeB1, dataB1['lum'].iloc[wc], mgeB2, dataB2['lum'].iloc[wc], tau)
                      i_logTB = inter_rout(mgeB1, dataB1['logT'].iloc[wc], mgeB2, dataB2['logT'].iloc[wc], tau)
                      i_MvB = inter_rout(mgeB1, dataB1['Mv'].iloc[wc], mgeB2, dataB2['Mv'].iloc[wc], tau)
                      i_jkB = inter_rout(mgeB1, dataB1['jk'].iloc[wc], mgeB2, dataB2['jk'].iloc[wc], tau)

                      i_min = inter_rout(bmet[j], i_minA, bmet[j+1], i_minB, metal)
                      i_mfin = inter_rout(bmet[j], i_mfinA, bmet[j+1], i_mfinB, metal)
                      i_lum = inter_rout(bmet[j], i_lumA, bmet[j+1], i_lumB, metal)
                      i_logT = inter_rout(bmet[j], i_logTA, bmet[j+1], i_logTB, metal)
                      i_Mv = inter_rout(bmet[j], i_MvA, bmet[j+1], i_MvB, metal)
                      i_jk = inter_rout(bmet[j], i_jkA, bmet[j+1], i_jkB, metal)
                      

                      mini[wc] = i_min
                      mfin[wc] = i_mfin
                      lum[wc] = i_lum
                      logT[wc] = i_logT
                      Mv[wc] = i_Mv
                      jk[wc] = i_jk
	
                   if order == 'metage':
                      in_min1 = inter_rout(bmet[j], dataA1['min'].iloc[wc], bmet[j+1], dataB1['min'].iloc[wc], metal)
                      in_fin1 = inter_rout(bmet[j], dataA1['mfin'].iloc[wc], bmet[j+1], dataB1['mfin'].iloc[wc], metal)
                      in_lum1 = inter_rout(bmet[j], dataA1['lum'].iloc[wc], bmet[j+1], dataB1['lum'].iloc[wc], metal)
                      in_logT1 = inter_rout(bmet[j], dataA1['logT'].iloc[wc], bmet[j+1], dataB1['logT'].iloc[wc], metal)
                      in_Mv1 = inter_rout(bmet[j], dataA1['Mv'].iloc[wc], bmet[j+1], dataB1['Mv'].iloc[wc], metal)
                      in_jk1 = inter_rout(bmet[j], dataA1['jk'].iloc[wc], bmet[j+1], dataB1['jk'].iloc[wc], metal)
		
                      in_min2 = inter_rout(bmet[j], dataA2['min'].iloc[wc], bmet[j+1], dataB2['min'].iloc[wc], metal)
                      in_fin2 = inter_rout(bmet[j], dataA2['mfin'].iloc[wc], bmet[j+1], dataB2['mfin'].iloc[wc], metal)
                      in_lum2 = inter_rout(bmet[j], dataA2['lum'].iloc[wc], bmet[j+1], dataB2['lum'].iloc[wc], metal)
                      in_logT2 = inter_rout(bmet[j], dataA2['logT'].iloc[wc], bmet[j+1], dataB2['logT'].iloc[wc], metal)
                      in_Mv2 = inter_rout(bmet[j], dataA2['Mv'].iloc[wc], bmet[j+1], dataB2['Mv'].iloc[wc], metal)
                      in_jk2 = inter_rout(bmet[j], dataA2['jk'].iloc[wc], bmet[j+1], dataB2['jk'].iloc[wc], metal)
			
                      in_min = inter_rout(mgeA1, in_min1, mgeA2, in_min2, tau)
                      in_fin = inter_rout(mgeA1, in_fin1, mgeA2, in_fin2, tau)
                      in_lum = inter_rout(mgeA1, in_lum1, mgeA2, in_lum2, tau)
                      in_logT = inter_rout(mgeA1, in_logT1, mgeA2, in_logT2, tau)
                      in_Mv = inter_rout(mgeA1, in_Mv1, mgeA2, in_Mv2, tau)
                      in_jk = inter_rout(mgeA1, in_jk1, mgeA2, in_jk2, tau)
	
                      mini[wc] = in_min
                      mfin[wc] = in_fin
                      lum[wc] = in_lum
                      logT[wc] = in_logT
                      Mv[wc] = in_Mv
                      jk[wc] = in_jk
	    					
            end = time.time()	
            #print(end-start)	            
            if Phot_sys == 'GAIA_EDR3':
               def read_iso_file(file_path):
	           
                   chunk = pd.read_csv(file_path, sep='\s+', comment='#', header=None,
                                      names=['min', 'mfin', 'lum', 'logT', 'MG', 'G_BP', 'G_RP'],
                                      usecols=[0, 1, 2, 3, 4, 5, 6],chunksize=200)
                   data = pd.concat(chunk)
                   return data 
               MODA1 = os.path.join(directory1, modA1)
               MODA2 = os.path.join(directory1, modA2)
               MODB1 = os.path.join(directory2, modB1)
               MODB2 = os.path.join(directory2, modB2)
	       
	
               dataA1 = read_iso_file(MODA1)
               dataA2 = read_iso_file(MODA2)
               dataB1 = read_iso_file(MODB1)
               dataB2 = read_iso_file(MODB2)

	       # Calculate G_GRP for each dataset
               dataA1['G_GRP'] = dataA1['G_BP'] - dataA1['G_RP']
               dataA2['G_GRP'] = dataA2['G_BP'] - dataA2['G_RP']
               dataB1['G_GRP'] = dataB1['G_BP'] - dataB1['G_RP']
               dataB2['G_GRP'] = dataB2['G_RP'] - dataB2['G_RP']	
	       
	       # Check for matching errors in the number of elements
               if len(dataA1) != len(dataA2):
                  print('MATCHING ERROR: dataA1 and dataA2 have different lengths')
               if len(dataB1) != len(dataB2):
                  print('MATCHING ERROR: dataB1 and dataB2 have different lengths')
               if len(dataA1) != len(dataB2):
                  print('MATCHING ERROR: dataA1 and dataB2 have different lengths')
	
	       # Define the number of elements in the isochrone
               nbasti = len(dataA1)
               
	       # Initialize arrays for storing interpolated values
               mini = np.zeros(nbasti)
               mfin = np.zeros(nbasti)
               lum = np.zeros(nbasti)
               logT = np.zeros(nbasti)
               MG = np.zeros(nbasti)
               G_GRP = np.zeros(nbasti)
	 
               for wc in range(nbasti):
                   if order == 'agemet':
                      i_minA = inter_rout(mgeA1, dataA1['min'].iloc[wc], mgeA2, dataA2['min'].iloc[wc], tau)
                      i_mfinA = inter_rout(mgeA1, dataA1['mfin'].iloc[wc], mgeA2, dataA2['mfin'].iloc[wc], tau)
                      i_lumA = inter_rout(mgeA1, dataA1['lum'].iloc[wc], mgeA2, dataA2['lum'].iloc[wc], tau)
                      i_logTA = inter_rout(mgeA1, dataA1['logT'].iloc[wc], mgeA2, dataA2['logT'].iloc[wc], tau)
                      i_MGA = inter_rout(mgeA1, dataA1['MG'].iloc[wc], mgeA2, dataA2['MG'].iloc[wc], tau)
                      i_G_GRPA = inter_rout(mgeA1, dataA1['G_GRP'].iloc[wc], mgeA2, dataA2['G_GRP'].iloc[wc], tau)
	
                      i_minB = inter_rout(mgeB1, dataB1['min'].iloc[wc], mgeB2, dataB2['min'].iloc[wc], tau)
                      i_mfinB = inter_rout(mgeB1, dataB1['mfin'].iloc[wc], mgeB2, dataB2['mfin'].iloc[wc], tau)
                      i_lumB = inter_rout(mgeB1, dataB1['lum'].iloc[wc], mgeB2, dataB2['lum'].iloc[wc], tau)
                      i_logTB = inter_rout(mgeB1, dataB1['logT'].iloc[wc], mgeB2, dataB2['logT'].iloc[wc], tau)
                      i_MGB = inter_rout(mgeB1, dataB1['MG'].iloc[wc], mgeB2, dataB2['MG'].iloc[wc], tau)
                      i_G_GRPB = inter_rout(mgeB1, dataB1['G_GRP'].iloc[wc], mgeB2, dataB2['G_GRP'].iloc[wc], tau)
	
                      i_min = inter_rout(bmet[j], i_minA, bmet[j+1], i_minB, metal)
                      i_mfin = inter_rout(bmet[j], i_mfinA, bmet[j+1], i_mfinB, metal)
                      i_lum = inter_rout(bmet[j], i_lumA, bmet[j+1], i_lumB, metal)
                      i_logT = inter_rout(bmet[j], i_logTA, bmet[j+1], i_logTB, metal)
                      i_MG = inter_rout(bmet[j], i_MGA, bmet[j+1], i_MGB, metal)
                      i_G_GRP = inter_rout(bmet[j], i_G_GRPA, bmet[j+1], i_G_GRPB, metal)
		
                      mini[wc] = i_min
                      mfin[wc] = i_mfin
                      lum[wc] = i_lum
                      logT[wc] = i_logT
                      MG[wc] = i_MG
                      G_GRP[wc] = i_G_GRP
	
                   if order == 'metage':
                      in_min1 = inter_rout(bmet[j], dataA1['min'].iloc[wc], bmet[j+1], dataB1['min'].iloc[wc], metal)
                      in_fin1 = inter_rout(bmet[j], dataA1['mfin'].iloc[wc], bmet[j+1], dataB1['mfin'].iloc[wc], metal)
                      in_lum1 = inter_rout(bmet[j], dataA1['lum'].iloc[wc], bmet[j+1], dataB1['lum'].iloc[wc], metal)
                      in_logT1 = inter_rout(bmet[j], dataA1['logT'].iloc[wc], bmet[j+1], dataB1['logT'].iloc[wc], metal)
                      in_MG1 = inter_rout(bmet[j], dataA1['MG'].iloc[wc], bmet[j+1], dataB1['MG'].iloc[wc], metal)
                      in_G_GRP1 = inter_rout(bmet[j], dataA1['G_GRP'].iloc[wc], bmet[j+1], dataB1['G_GRP'].iloc[wc], metal)
		
                      in_min2 = inter_rout(bmet[j], dataA2['min'].iloc[wc], bmet[j+1], dataB2['min'].iloc[wc], metal)
                      in_fin2 = inter_rout(bmet[j], dataA2['mfin'].iloc[wc], bmet[j+1], dataB2['mfin'].iloc[wc], metal)
                      in_lum2 = inter_rout(bmet[j], dataA2['lum'].iloc[wc], bmet[j+1], dataB2['lum'].iloc[wc], metal)
                      in_logT2 = inter_rout(bmet[j], dataA2['logT'].iloc[wc], bmet[j+1], dataB2['logT'].iloc[wc], metal)
                      in_MG2 = inter_rout(bmet[j], dataA2['MG'].iloc[wc], bmet[j+1], dataB2['MG'].iloc[wc], metal)
                      in_G_GRP2 = inter_rout(bmet[j], dataA2['G_GRP'].iloc[wc], bmet[j+1], dataB2['G_GRP'].iloc[wc], metal)
			
                      in_min = inter_rout(mgeA1, in_min1, mgeA2, in_min2, tau)
                      in_fin = inter_rout(mgeA1, in_fin1, mgeA2, in_fin2, tau)
                      in_lum = inter_rout(mgeA1, in_lum1, mgeA2, in_lum2, tau)
                      in_logT = inter_rout(mgeA1, in_logT1, mgeA2, in_logT2, tau)
                      in_MG = inter_rout(mgeA1, in_MG1, mgeA2, in_MG2, tau)
                      in_G_GRP = inter_rout(mgeA1, in_G_GRP1, mgeA2, in_G_GRP2, tau)
	
                      mini[wc] = in_min
                      mfin[wc] = in_fin
                      lum[wc] = in_lum
                      logT[wc] = in_logT
                      MG[wc] = in_MG
                      G_GRP[wc] = in_G_GRP
	    
            #print(dataA1['min'])
            #print(mgeA1,mgeA2,tau)
	    # Compute log(g) at each point along the isochrone
            lumS = 10.**lum
            Tmod = 10.**logT
            R2 = (Tsun/Tmod)**4 * lumS
 
            logg = logsun + np.log10(mini/R2)  # log(g)
            
 
            min_mass = np.min(mini)
            max_mass = np.max(mini)
			
            
            def randomp(pow, n, range_x=[5, 100], seed=None):
                """
	    	Generates an array of random numbers distributed as a power law.
	    
	    	Parameters:
	                pow (float): Exponent of power law.
	                n (int): Number of elements in the generated vector.
	                range_x (list or tuple, optional): 2-element vector [low, high] specifying the range of 
	                output X valuees. Default is [5, 100].
	                seed (int, optional): Seed value for numpy.random.default_rng. Default is None.
                PROCEDURE:  
                        "Transformation Method" for random variables is described in Bevington 
                         & Robinson, "Data Reduction & Error Analysis for Physical Sciences", 2nd
                         Edition (McGraw-Hill, 1992). p. 83.
	    	Returns:
	                x (ndarray): Vector of random numbers, distributed as a power law between specified range.
                """
	    	# Set seed if provided
                if seed is not None:
                   np.random.seed(seed)
	    
	    	# Unpack range
                lo, hi = sorted(range_x)
	    
	    	# Generate random numbers
                r = np.random.rand(n)
	        
                pow1 = pow + 1
	    	# Compute random numbers according to the power law
                if pow != -1.0:
                   norm = 1.0 / (hi**pow1 - lo**pow1)
                   expo = np.log10(r / norm + lo**pow1) / pow1
                   x = 10.0 ** expo
                else:
                   norm = 1.0 / (np.log(hi) - np.log(lo))
                   x = np.exp(r / norm + np.log(lo))
	    
                return x
			
	    # Populate isochrone with IMF
	    # Assuming randomp is a function that returns a random sample from an IMF
            #print(min_mass,max_mass)
            if self.imf_type == 'salpeter':
               mass = randomp(self.imf_slope, self.nstars, range_x=[min_mass, max_mass])
               Nstars = len(mass)
	    
            if self.imf_type == 'kroupa' or self.imf_type == 'chabrierlognormal':
               mass = imf.inverse_imf(np.random.random(self.nstars), min_mass, max_mass, massfunc= self.imf_type) # 'chabrierlognormal' or 'kroupa'
               Nstars = len(mass)
	        
            if Phot_sys == '2mass':
	    
               # Interpolate Mk, jk, and logg to the generated masses
               Mv_interpol = interp1d(mini, Mv, kind='linear', fill_value="extrapolate")(mass)
               jk_interpol = interp1d(mini, jk, kind='linear', fill_value="extrapolate")(mass)
               lg_interpol = interp1d(mini, logg, kind='linear', fill_value="extrapolate")(mass)
	            
               color_idx = jk_interpol
               abs_mag = Mv_interpol
               iso_logg = lg_interpol
               
            if Phot_sys == 'GAIA_EDR3':
	    
               # Interpolate MG, G_GRP, and logg to the generated masses
               MG_interpol = interp1d(mini, MG, kind='linear', fill_value="extrapolate")(mass)
               G_GRP_interpol = interp1d(mini, G_GRP, kind='linear', fill_value="extrapolate")(mass)
               lg_interpol = interp1d(mini, logg, kind='linear', fill_value="extrapolate")(mass)
	            
               color_idx = G_GRP_interpol
               abs_mag = MG_interpol
               iso_logg = lg_interpol   
                       	
            return color_idx, abs_mag, iso_logg	,Nstars									
