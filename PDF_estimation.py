import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from scipy import stats, interpolate
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import dblquad
from scipy.integrate import quad
import math
from astropy.io import fits
from scipy import stats
import sys

# Access astronomical databases
from pyvo import registry  # version >=1.4.1 

# Moc and HEALPix tools
from mocpy import MOC

# Sky visualization
from ipyaladin import Aladin    # version >=0.3.0


#**(270,-73)
#(270,4)
#***(69,-38)
#(67,-13)
#(4.5,-72.5)
#(63,-6)
#*(12,-63)

# specific_feh = -0.5  # Example value for feh -1.5
# specific_distance = 5200  # Example value for distance 210
long = float(sys.argv[1])                                 
lat = float(sys.argv[2])
conesearch_radius = 1.78
path = 'galaxia.dat'
path_iso = 'test_results.dat'
bin_size_distance = 200
bin_size_feh = 50
bin_size_age = 40
bin_along_distance = False 
bin_along_distance_feh = True if sys.argv[3] == 'True' else False
bin_along_distance_age = True if sys.argv[4] == 'True' else False


# the catalogue name in VizieR
CATALOGUE = "J/MNRAS/478/4513"

# each resource in the VO has an identifier, called ivoid. For vizier catalogs,
# the VO ids can be constructed like this:
catalogue_ivoid = f"ivo://CDS.VizieR/{CATALOGUE}"
# the actual query to the registry
voresource = registry.search(ivoid=catalogue_ivoid)[0]

conesearch_center = (long, lat)
conesearch_records = voresource.get_service("conesearch").search(
    pos=conesearch_center,
    sr=conesearch_radius,)

Star_ID = conesearch_records['StarId'].data[0:len(conesearch_records)]

hdul = fits.open('GALAH_DR3_main_allstar_v2.fits')
data_obs = hdul[1].data

# Create a dictionary for fast lookup of indices in data_obs
obs_dict = {row[0]: idx for idx, row in enumerate(data_obs)}

# Use the dictionary to find indices for Star_ID
index = [obs_dict[star] for star in Star_ID if star in obs_dict]


data_obs = pd.DataFrame(data_obs)

index_flagged = []
for i in index:
    if data_obs['flag_sp'][i] == 0 and data_obs['flag_fe_h'][i] == 0:
        index_flagged.append(i)

feh = []
logg = []
teff = []
ebv = [] 
e_fe_h = []
sobject_id = []
for i in index_flagged:
    feh.append(data_obs['fe_h'][i])
    logg.append(data_obs['logg'][i])
    teff.append(data_obs['teff'][i])
    ebv.append(data_obs['ebv'][i])
    sobject_id.append(data_obs['sobject_id'][i])
    e_fe_h.append(data_obs['e_fe_h'])

# print(sobject_id)

hdul = fits.open('GALAH_DR3_VAC_ages_v2.fits')

data_gaiadr3 = pd.DataFrame(hdul[1].data)
# Ensure all columns in the dataframe are little-endian
# Apply little-endian conversion only to numeric columns (int and float)
data_gaiadr3 = data_gaiadr3.apply(lambda col: col.astype('<f8') if col.dtype.kind == 'f' else 
                                   (col.astype('<i8') if col.dtype.kind == 'i' else col))
data_gaiadr3 = data_gaiadr3.dropna()

sobject_id = np.array(sobject_id, dtype=np.int64)
# print(sobject_id.dtype)  # Check the dtype of the array
# print(data_gaiadr3['sobject_id'].dtype)  # Check the dtype of the column
index_gaiadr3 = data_gaiadr3[data_gaiadr3['sobject_id'].isin(sobject_id)].index.tolist()


age = []
distance = []
ebv_bstep = []
e_age_bstep = []
e50_age_bstep = []
for i in index_gaiadr3:
    age.append(data_gaiadr3['age_bstep'][i])
    distance.append(1/((data_gaiadr3['distance_bstep'][i])*pow(10,-3)))
    ebv_bstep.append(data_gaiadr3['ebv_bstep'][i])
    e_age_bstep.append(data_gaiadr3['e_age_bstep'][i])
    e50_age_bstep.append(data_gaiadr3['e50_age_bstep'][i])


filtered_feh = [f for i,f in enumerate(feh) if -3.26 <= f <= 0.5 and 0.5 <= age[i] <= 14.375]
filtered_age = [a for i,a in enumerate(age) if 0.5 <= a <= 14.375 and -3.26 <= feh[i] <= 0.5]
filtered_distance = [d for i,d in enumerate(distance) if -3.26 <= feh[i] <= 0.5]

feh_distance = list(zip(filtered_feh, filtered_distance))
age_distance = list(zip(filtered_age, filtered_distance))

#print(feh_distance)
#print(age_distance)

# # Check for NaN values before conversion
# print(data_gaiadr3.isna().sum())

# # Check the unique values in the columns causing issues
# print(data_gaiadr3['sobject_id'].unique())

# # Check for non-numeric values
# for col in data_gaiadr3.columns:
#     if data_gaiadr3[col].dtype not in ['float64', 'int64']:
#         print(f"Column {col} contains non-numeric data.")

# # Check the unique values for a few columns that have NaN
# print(data_gaiadr3['age_bstep'].unique())
# # Checking for any invalid (or missing) entries in the relevant columns
# print(data_gaiadr3['age_bstep'].isna().sum())  # Count NaN values in this column
# print(data_gaiadr3['age_bstep'].head())  # Show a preview of the column


age_peak = []
feh_peak = []
age_median_peak = []
feh_median_peak = []
age_16_percent = []
age_84_percent = []
feh_16_percent = []
feh_84_percent = []

def read_columns(filepath):
    
    chunk = pd.read_csv(filepath, sep = '\s+', comment='#', header=None,
                         names = ['mini','mfin','age','feh','distance','l','b','px','py','pz','popid','exbv_schlegel','log_g'],
                         usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12], chunksize = 200)
    data = pd.concat(chunk)
    
    return data
    
def calculate_conditional_pdf_age(parameter_values):
    # age_pdf_values = age_kde(age_values)
    
    def integrate_kde(kde, lower, upper):
        result, _ = quad(kde, lower, upper)
        return result
    parameter_interval = parameter_values[1]-parameter_values[0]
    probability = []

    for i in range(len(parameter_values)):
        probability.append(integrate_kde(age_kde,parameter_values[i]-parameter_interval,parameter_values[i]))
    return np.array(probability)  
    
def calculate_conditional_pdf_feh(parameter_values):
    # age_pdf_values = age_kde(age_values)
    
    def integrate_kde(kde, lower, upper):
        result, _ = quad(kde, lower, upper)
        return result
    parameter_interval = parameter_values[1]-parameter_values[0]
    probability = []

    for i in range(len(parameter_values)):
        probability.append(integrate_kde(feh_kde,parameter_values[i]-parameter_interval,parameter_values[i]))
    return np.array(probability)   
    
def calculate_conditional_pdf(parameter_values):
    # age_pdf_values = age_kde(age_values)
    
    def integrate_kde(kde, lower, upper):
        result, _ = quad(kde, lower, upper)
        return result
    parameter_interval = parameter_values[1]-parameter_values[0]
    probability = []

    if bin_along_distance_feh == True:
        for i in range(len(parameter_values)):
            probability.append(integrate_kde(age_kde,parameter_values[i]-parameter_interval,parameter_values[i]))
    if bin_along_distance_age == True:
        for i in range(len(parameter_values)):
            probability.append(integrate_kde(feh_kde,parameter_values[i]-parameter_interval,parameter_values[i]))        
    return np.array(probability)

def read_iso_file(file_path):
	# Read the file into a pandas DataFrame
    chunk = pd.read_csv(file_path, sep='\s+', comment='#', header=None,
                        names=['Age', 'feh', 'distance', 'prob'],
                        usecols=[0, 1, 2, 3],chunksize=200)
    data_iso = pd.concat(chunk)
    return data_iso

def inter_rout(x1, y1, x2, y2, x):
    y = ( (y2-y1)/(x2-x1) ) * (x-x1) + y1
    return y
    
data = read_columns(path)
data_iso = read_iso_file(path_iso)

count_age_11 = 0
count_age_14 = 0
count_age_10 = 0
Age = []
index_model = []
Age_model = []
feh_model = []


for i in range(len(data['age'])): 
    index_model.append(i)

for index in index_model:
    Age_model.append(data['age'][index])

for i in range(len(Age_model)):
    
    if Age_model[i] == 11.000:
        count_age_11+=1
    elif Age_model[i] == 14.000:
        count_age_14+=1
    elif Age_model[i] ==10.000 and (Age_model[i-1] == 14.000 or Age_model[i-1] == 10.000):
        count_age_10+=1
    else:
        Age.append(Age_model[i])
 
age_11 = stats.norm.rvs(loc=11, scale =0.5, size = count_age_11)
for i in range(len(age_11)):
    Age.append(age_11[i])

if count_age_10 != 0:
    
    age_14 = stats.norm.rvs(loc=14, scale =0.5, size = count_age_14*2)
    count=0
    for i in range(len(age_14)):
        age_14.sort()
        count+=1
        if count <= count_age_14:
            Age.append(age_14[i])
    
    age_10 = stats.norm.rvs(loc=10, scale =0.5, size = count_age_10)
    for i in range(len(age_10)):
        Age.append(age_10[i]) 
else:
    
    age_14 = stats.norm.rvs(loc=14, scale =0.5, size = count_age_14*2)
    count=0
    for i in range(len(age_14)):
        age_14.sort()
        count+=1
        if count <= count_age_14:
            Age.append(age_14[i]) 
        

total_stars = len(data)
stars = {}
nstar_dist_bin = int(total_stars/bin_size_distance)
data['age'] = np.array(Age)

data.sort_values(by=['distance'], ascending=[True], inplace=True, ignore_index=True)

count = 1
for i in range(total_stars):
    
    try:
        value = stars["stars_distance_"+str(count)]
    except KeyError:
        stars["stars_distance_"+str(count)] = []

    if (count-1)*bin_size_distance <= len(stars["stars_distance_" + str(count)]) + (count-1)*bin_size_distance <= count*bin_size_distance:
        if bin_along_distance_feh == True or bin_along_distance == True:
            (stars["stars_distance_"+str(count)]).append([data['feh'][i],data['distance'][i],data['age'][i]])
        if bin_along_distance_age == True:
            (stars["stars_distance_"+str(count)]).append([data['age'][i],data['distance'][i],data['feh'][i]])
            
    if len(stars["stars_distance_" + str(count)]) == bin_size_distance:
        stars["stars_distance_"+str(count)].sort()
        count+=1

if bin_along_distance_feh == True:
    print('Heyy :)')
    stars_fully_binned = {}

    for count in range(1,len(stars)+1):
        count_feh = 1
        for i in range(len(stars["stars_distance_"+str(count)])):
            try:
                value = stars_fully_binned["stars_distance_"+str(count)+'_feh_'+str(count_feh)]
            except KeyError:
                stars_fully_binned["stars_distance_"+str(count)+'_feh_'+str(count_feh)] = []

            if (count_feh-1)*bin_size_feh <= len(stars_fully_binned["stars_distance_"+str(count)+'_feh_'+str(count_feh)]) + (count_feh-1)*bin_size_feh <= count_feh*bin_size_feh:
                (stars_fully_binned["stars_distance_"+str(count)+'_feh_'+str(count_feh)]).append(tuple(stars.items())[count-1][1][i])
        
            if len(stars_fully_binned["stars_distance_"+str(count)+'_feh_'+str(count_feh)]) == bin_size_feh:
                count_feh+=1

if bin_along_distance_age == True:

    stars_fully_binned = {}

    for count in range(1,len(stars)+1):
        count_age = 1
        for i in range(len(stars["stars_distance_"+str(count)])):
            try:
                value = stars_fully_binned["stars_distance_"+str(count)+'_age_'+str(count_age)]
            except KeyError:
                stars_fully_binned["stars_distance_"+str(count)+'_age_'+str(count_age)] = []

            if (count_age-1)*bin_size_age <= len(stars_fully_binned["stars_distance_"+str(count)+'_age_'+str(count_age)]) + (count_age-1)*bin_size_age <= count_age*bin_size_age:
                (stars_fully_binned["stars_distance_"+str(count)+'_age_'+str(count_age)]).append(tuple(stars.items())[count-1][1][i])
        
            if len(stars_fully_binned["stars_distance_"+str(count)+'_age_'+str(count_age)]) == bin_size_age:
                count_age+=1
                
survey_prob = []
survey_prob_model = []
if bin_along_distance == True:
    for (specific_feh,specific_distance) in feh_distance:

        for i in range(len(data)):
            if specific_distance < data['distance'][i]:
                count_distance = int(i/bin_size_distance) + 1
                break
        Age_data = []
        feh_data = []
        distance_data = []

        for i in range(len(stars['stars_distance_'+str(count_distance)])):
            Age_data.append(stars['stars_distance_'+str(count_distance)][i][2]) 
            distance_data.append(stars['stars_distance_'+str(count_distance)][i][1]) 
            feh_data.append(stars['stars_distance_'+str(count_distance)][i][0])
        
        data_iso_prob = data_iso['prob'] 
        data_iso_age = data_iso['Age']
        data_iso_feh = data_iso['feh']
        data_iso_distance = data_iso['distance']
        age_kde = gaussian_kde(Age_data)        
        feh_kde = gaussian_kde(feh_data)

        distance_count = 1
        for i in range(1,len(data_iso_feh)):
            if data_iso_feh[i] == data_iso_feh[i-1]:
                distance_count+=1
            else:
                break
        age_1 = []
        prob_1=[]
        for i in range(len(data_iso_feh)):

            if data_iso_feh[i] == specific_feh <= data_iso_feh[i+1] and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:
            
                age_1.append(inter_rout(data_iso_distance[i],data_iso_age[i],data_iso_distance[i+1],data_iso_age[i+1],specific_distance))
                prob_1.append(inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance))
        
    
            elif data_iso_feh[i] < specific_feh < data_iso_feh[i+distance_count] and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:

                age1 = inter_rout(data_iso_distance[i],data_iso_age[i],data_iso_distance[i+1],data_iso_age[i+1],specific_distance)
                prob1 = inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance)

                age2 = inter_rout(data_iso_distance[i+distance_count],data_iso_age[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_age[i+distance_count+1],specific_distance)
                prob2 = inter_rout(data_iso_distance[i+distance_count],data_iso_prob[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_prob[i+distance_count+1],specific_distance)

                age_1.append(inter_rout(data_iso_feh[i],age1,data_iso_feh[i+distance_count],age2,specific_feh))
                prob_1.append(inter_rout(data_iso_feh[i],prob1,data_iso_feh[i+distance_count],prob2,specific_feh))

        age_1 = np.array(age_1)

        distance_count = 1
        for i in range(1,len(data_iso_age)):
            if data_iso_age[i] == data_iso_age[i-1]:
                distance_count+=1
            else:
                break

        feh_1 = []
        prob_2=[]
        for i in range(len(data_iso_age)):
            if data_iso_age[i] == specific_age and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:

                feh_1.append(inter_rout(data_iso_distance[i],data_iso_feh[i],data_iso_distance[i+1],data_iso_feh[i+1],specific_distance))
                prob_2.append(inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance))

            elif data_iso_age[i] < specific_age < data_iso_age[i+distance_count] and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:

                feh1 = inter_rout(data_iso_distance[i],data_iso_feh[i],data_iso_distance[i+1],data_iso_feh[i+1],specific_distance)
                prob1 = inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance)

                feh2 = inter_rout(data_iso_distance[i+distance_count],data_iso_feh[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_feh[i+distance_count+1],specific_distance)
                prob2 = inter_rout(data_iso_distance[i+distance_count],data_iso_prob[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_prob[i+distance_count+1],specific_distance)

                feh_1.append(inter_rout(data_iso_age[i],feh1,data_iso_age[i+distance_count],feh2,specific_age))
                prob_2.append(inter_rout(data_iso_age[i],prob1,data_iso_age[i+distance_count],prob2,specific_age))
        

        #if specific_feh == daTrueta_iso_feh[i] and specific_distance == data_iso_distance[i]:
        #   index.append(i)
        feh_1 = np.array(feh_1) 
        #plt.plot(age,conditional_pdf_values,label='Model probablilty Age')
        #plt.plot(age,conditional_pdf_values*prob,label='Observed Selection Age')
        #plt.plot(age,prob,label="Selection_prob")                                                                                


        #plt.xlabel('Age')
        #plt.ylabel('Probability of Observing a Star')
        #plt.legend()
        #plt.grid()
        #plt.show()


survey_prob = []
if bin_along_distance_feh == True:
    count_age = 0
    for (specific_feh,specific_distance) in feh_distance:

        for i in range(len(data)):
            if specific_distance < data['distance'][i]:
                count_distance = int(i/bin_size_distance) + 1
                break

        if bin_size_distance%bin_size_feh == 0:
            bin_ratio = int(bin_size_distance/bin_size_feh)
        else:
            bin_ratio = int(bin_size_distance/bin_size_feh)+1
    
        for Bin in range(1,bin_ratio,1):
            for i in range(len(stars_fully_binned['stars_distance_'+str(count_distance)+'_feh_'+str(Bin)])):
                if stars_fully_binned['stars_distance_'+str(count_distance)+'_feh_'+str(Bin)][i][0] > specific_feh:
                    count_feh = Bin
                    break
            else:
                continue
            break
        Age_data = []
        feh_data = []
        distance_data = []

        for i in range(len(stars_fully_binned['stars_distance_'+str(count_distance)+'_feh_'+str(count_feh)])):
            Age_data.append(stars_fully_binned['stars_distance_'+str(count_distance)+'_feh_'+str(count_feh)][i][2]) 
            distance_data.append(stars_fully_binned['stars_distance_'+str(count_distance)+'_feh_'+str(count_feh)][i][1]) 
            feh_data.append(stars_fully_binned['stars_distance_'+str(count_distance)+'_feh_'+str(count_feh)][i][0]) 

        data_iso_prob = data_iso['prob'] 
        data_iso_age = data_iso['Age']
        data_iso_feh = data_iso['feh']
        data_iso_distance = data_iso['distance']
        age_kde = gaussian_kde(Age_data)

        distance_count = 1
        for i in range(1,len(data_iso_feh)):
            if data_iso_feh[i] == data_iso_feh[i-1]:
                distance_count+=1
            else:
                break
        age = []
        prob=[]
        for i in range(len(data_iso_feh)):

            if data_iso_feh[i] == specific_feh <= data_iso_feh[i+1] and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:
            
                age.append(inter_rout(data_iso_distance[i],data_iso_age[i],data_iso_distance[i+1],data_iso_age[i+1],specific_distance))
                prob.append(inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance))
        
    
            elif data_iso_feh[i] < specific_feh < data_iso_feh[i+distance_count] and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:

                age1 = inter_rout(data_iso_distance[i],data_iso_age[i],data_iso_distance[i+1],data_iso_age[i+1],specific_distance)
                prob1 = inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance)

                age2 = inter_rout(data_iso_distance[i+distance_count],data_iso_age[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_age[i+distance_count+1],specific_distance)
                prob2 = inter_rout(data_iso_distance[i+distance_count],data_iso_prob[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_prob[i+distance_count+1],specific_distance)

                age.append(inter_rout(data_iso_feh[i],age1,data_iso_feh[i+distance_count],age2,specific_feh))
                prob.append(inter_rout(data_iso_feh[i],prob1,data_iso_feh[i+distance_count],prob2,specific_feh))
    
        age = np.array(age)   
        prob = np.array(prob)
        conditional_pdf_values = calculate_conditional_pdf(age)
        survey_prob.append(conditional_pdf_values*prob)
        survey_prob_model.append(conditional_pdf_values)
        conditional_pdf_values_prob = conditional_pdf_values*prob
        for i in range(len(conditional_pdf_values_prob)):
            if conditional_pdf_values_prob[i] == np.max(conditional_pdf_values_prob):
                break
        age_peak.append(age[i])
        sum_pdf = np.sum(conditional_pdf_values_prob)
        
        temp = 0
        for i in range(len(conditional_pdf_values_prob)):
            temp = temp + conditional_pdf_values_prob[i]
            if temp > sum_pdf/2:
                break
        age_median_peak.append(age[i])

        temp = 0
        for i in range(len(conditional_pdf_values_prob)):
            temp = temp + conditional_pdf_values_prob[i]
            if temp > sum_pdf/6.25:
                break
        age_16_percent.append(age[i])

        temp = 0
        for i in range(len(conditional_pdf_values_prob)):
            temp = temp + conditional_pdf_values_prob[i]
            if temp > sum_pdf/1.19:
                break
        age_84_percent.append(age[i])
        
        print(f"The Age PDF of star at ([Fe/H],Distance)=({specific_feh:.2f},{specific_distance:.1f}) has been estimated") 
    
    x = np.zeros(len(survey_prob[0]))
    for i in range(len(survey_prob)):
        x = x + survey_prob[i]
    
    final_prob = x/(np.trapz(x,age))
    #print(final_prob)
    # data_age = 

    y = np.zeros(len(survey_prob_model[0]))
    for i in range(len(survey_prob_model)):
        y = y + survey_prob_model[i]
    final_prob_model = y/(np.trapz(y,age))
    
    # # 1. Calculate the ECDF for the observed data
    # data_age = [a for a,d in age_distance]
    # data_sorted = np.sort(data_age)
    # ecdf_data = np.arange(1, len(data_sorted) + 1) / len(data_sorted)

    # # 2. Calculate the theoretical CDF
    # cdf_theoretical = np.cumsum(final_prob)  # Cumulative sum of theoretical probabilities

    # # 3. Interpolate the theoretical CDF at the observed data points
    # # This step is to align the observed data with the theoretical CDF
    # cdf_interpolated_theoretical = np.interp(data_sorted, age, cdf_theoretical)

    # # 4. Perform the KS test (comparing ECDF of data and theoretical CDF)
    # ks_statistic, p_value = stats.ks_2samp(ecdf_data, cdf_interpolated_theoretical)

    # # Output the results
    # print(f"KS Statistic: {ks_statistic}")
    # print(f"P-Value: {p_value}")

    # Sort the data and theoretical distributions
    theoretical_cdf = np.cumsum(final_prob)
    theoretical_cdf /= theoretical_cdf[-1]  # Normalize to ensure the total is 1
    
    # Generate quantiles for the theoretical distribution
    theoretical_quantiles = np.interp(
        np.linspace(0, 1, len(age)), theoretical_cdf, age
    )
    
    data_age = [a for a,d in age_distance]
    
    # Sort data ages to find quantiles
    data_quantiles = np.percentile(data_age, np.linspace(0, 100, len(age)))
    
    # Plot the Q-Q plot
    plt.figure(figsize=(8, 6),dpi=300)
    plt.scatter(theoretical_quantiles, data_quantiles, s=5, color='orange', label='Data vs. Framework Model')
    plt.plot(theoretical_quantiles, theoretical_quantiles, 'b--', label='Perfect Match (y=x)')
    plt.title('Q-Q Plot: Data vs. Framework Model Age Distribution')
    plt.xlabel('Framework Model Quantiles')
    plt.ylabel('Data Quantiles')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'/Users/advik/OPACOS/RESULTS/Data Age vs Framework Model Age Q-Q plot for ({long},{lat}).png',dpi=300)
    #plt.show()
        
    # print(age_peak,age_median_peak,age_16_percent,age_84_percent)
    plt.figure(figsize=(7,14))
    plt.subplot(211)
    plt.plot(age,final_prob,label='Framework Age', color= 'blue')
    plt.xlabel('Age')
    plt.ylabel('Probability of Observing a Star with our Framework')    
    plt.subplot(212)    
    plt.plot(age,final_prob_model,label='GALAXIA Age', color = 'red')
    plt.xlabel('Age')
    plt.ylabel('Probability of Spawning a Star using GALAXIA')
    plt.legend()
    plt.grid()
    # plt.scatter(age_peak,gaia_age)
    plt.savefig(f'/Users/advik/OPACOS/RESULTS/Framework Inferred Age and GALAXIA Age PDF for ({long},{lat}).png')
    #plt.show()

    data_age = pd.DataFrame(filtered_age)
    # random number generator

    # Plot pandas histogram from dataframe with df.plot.hist (not df.hist)
    ax = data_age.plot.hist(bins=30, density=True, edgecolor='red', linewidth=0.125)

    # Save default x-axis limits for final formatting because the pandas kde
    # plot uses much wider limits which usually decreases readability
    xlim = ax.get_xlim()

    # Plot pandas KDE
    data_age.plot.density(color='blue', alpha=1, ax=ax) # same as df['var'].plot.kde()

    # Reset x-axis limits and edit legend and add title
    ax.set_xlim(xlim)
    ax.legend(labels=['KDE'], frameon=False)
    ax.set_title('Histogram overlaid with KDE', fontsize=14, pad=15)
    plt.savefig(f'/Users/advik/OPACOS/RESULTS/Data Age Histogram for ({long},{lat}).png')

survey_prob = []
if bin_along_distance_age == True:
    for (specific_age,specific_distance) in age_distance:
        
        for i in range(len(data)):
            if specific_distance < data['distance'][i]:
                count_distance = int(i/bin_size_distance) + 1
                break

        if bin_size_distance%bin_size_age == 0:
            bin_ratio = int(bin_size_distance/bin_size_age)
        else:
            bin_ratio = int(bin_size_distance/bin_size_age)+1
            
        for Bin in range(1,bin_ratio,1):
            for i in range(len(stars_fully_binned['stars_distance_'+str(count_distance)+'_age_'+str(Bin)])):
                if stars_fully_binned['stars_distance_'+str(count_distance)+'_age_'+str(Bin)][i][0] > specific_age:
                    count_age = Bin
                    break
            else:
                continue
            break
        Age_data = []
        feh_data = []
        distance_data = []

        for i in range(len(stars_fully_binned['stars_distance_'+str(count_distance)+'_age_'+str(count_age)])):
            Age_data.append(stars_fully_binned['stars_distance_'+str(count_distance)+'_age_'+str(count_age)][i][0]) 
            distance_data.append(stars_fully_binned['stars_distance_'+str(count_distance)+'_age_'+str(count_age)][i][1]) 
            feh_data.append(stars_fully_binned['stars_distance_'+str(count_distance)+'_age_'+str(count_age)][i][2]) 

        data_iso_prob = data_iso['prob'] 
        data_iso_age = data_iso['Age']
        data_iso_feh = data_iso['feh']
        data_iso_distance = data_iso['distance']
        feh_kde = gaussian_kde(feh_data)

        distance_count = 1
        for i in range(1,len(data_iso_age)):
            if data_iso_age[i] == data_iso_age[i-1]:
                distance_count+=1
            else:
                break

        feh = []
        prob=[]
        for i in range(len(data_iso_age)):
            if data_iso_age[i] == specific_age and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:

                feh.append(inter_rout(data_iso_distance[i],data_iso_feh[i],data_iso_distance[i+1],data_iso_feh[i+1],specific_distance))
                prob.append(inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance))

            elif data_iso_age[i] < specific_age < data_iso_age[i+distance_count] and data_iso_distance[i] <= specific_distance < data_iso_distance[i+1]:

                feh1 = inter_rout(data_iso_distance[i],data_iso_feh[i],data_iso_distance[i+1],data_iso_feh[i+1],specific_distance)
                prob1 = inter_rout(data_iso_distance[i],data_iso_prob[i],data_iso_distance[i+1],data_iso_prob[i+1],specific_distance)

                feh2 = inter_rout(data_iso_distance[i+distance_count],data_iso_feh[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_feh[i+distance_count+1],specific_distance)
                prob2 = inter_rout(data_iso_distance[i+distance_count],data_iso_prob[i+distance_count],data_iso_distance[i+distance_count+1],data_iso_prob[i+distance_count+1],specific_distance)

                feh.append(inter_rout(data_iso_age[i],feh1,data_iso_age[i+distance_count],feh2,specific_age))
                prob.append(inter_rout(data_iso_age[i],prob1,data_iso_age[i+distance_count],prob2,specific_age))


        #if specific_feh == daTrueta_iso_feh[i] and specific_distance == data_iso_distance[i]:
        #   index.append(i)
        feh = np.array(feh)
        prob = np.array(prob)
        
        conditional_pdf_values = calculate_conditional_pdf(feh)        
        survey_prob.append(conditional_pdf_values*prob)
        survey_prob_model.append(conditional_pdf_values)
        conditional_pdf_values_prob = conditional_pdf_values*prob
        for i in range(len(conditional_pdf_values_prob)):
            if conditional_pdf_values_prob[i] == np.max(conditional_pdf_values_prob):
                break         
        feh_peak.append(feh[i])
        sum_pdf = np.sum(conditional_pdf_values_prob)

        temp = 0
        for i in range(len(conditional_pdf_values_prob)):
            temp = temp + conditional_pdf_values_prob[i]
            if temp > sum_pdf/2:
                break
        feh_median_peak.append(feh[i])

        temp = 0
        for i in range(len(conditional_pdf_values_prob)):
            temp = temp + conditional_pdf_values_prob[i]
            if temp > sum_pdf/6.25:
                break
        feh_16_percent.append(feh[i])

        temp = 0
        for i in range(len(conditional_pdf_values_prob)):
            temp = temp + conditional_pdf_values_prob[i]
            if temp > sum_pdf/1.19:
                break
        feh_84_percent.append(feh[i])
        #plt.plot(age,conditional_pdf_values,label='Model probablilty Age')
        #plt.plot(age,conditional_pdf_values*prob,label='Observed Selection Age')
        #plt.plot(age,prob,label="Selection_prob")

        #plt.xlabel('Age')
        #plt.ylabel('Probability of Observing a Star')
        #plt.legend()
        #plt.grid()
        #plt.show()
#x=np.sum(survey_prob[0]*survey_prob[1]*survey_prob[2]*survey_prob[3])
#final_prob = (survey_prob[0]*survey_prob[1]*survey_prob[2]*survey_prob[3])/x
        print(f"The [Fe/H] PDF of star at (Age,Distance)=({specific_age:.2f},{specific_distance:.1f}) has been estimated") 

    x = np.zeros(len(survey_prob[0]))
    for i in range(len(survey_prob)):
        x = x + survey_prob[i]
    final_prob = x/(np.trapz(x,feh))

    y = np.zeros(len(survey_prob_model[0]))
    for i in range(len(survey_prob_model)):
        y = y + survey_prob_model[i]
    final_prob_model = y/(np.trapz(y,feh))
    
    
    # Sort the data and theoretical distributions
    theoretical_cdf = np.cumsum(final_prob)
    theoretical_cdf /= theoretical_cdf[-1]  # Normalize to ensure the total is 1
    
    # Generate quantiles for the theoretical distribution
    theoretical_quantiles = np.interp(
        np.linspace(0, 1, len(feh)), theoretical_cdf, feh
    )
    
    data_feh = [f for f,d in feh_distance]
    
    # Sort data ages to find quantiles
    data_quantiles = np.percentile(data_feh, np.linspace(0, 100, len(feh)))
    
    # Plot the Q-Q plot
    plt.figure(figsize=(8, 6),dpi=300)
    plt.scatter(theoretical_quantiles, data_quantiles, s=5, color='orange', label='Data vs. Framework Model')
    plt.plot(theoretical_quantiles, theoretical_quantiles, 'b--', label='Perfect Match (y=x)')
    plt.title('Q-Q Plot: Data vs. Framework Model [Fe/H] Distribution')
    plt.xlabel('Framework Model Quantiles')
    plt.ylabel('Data Quantiles')
    plt.legend()
    plt.grid(True)
    fe_h = "Fe-H"
    plt.savefig(f'/Users/advik/OPACOS/RESULTS/Data {fe_h} vs Framework Model {fe_h} Q-Q plot for ({long},{lat}).png',dpi=300)
    #plt.show()
    
    #print(feh_peak,feh_median_peak,feh_16_percent,feh_84_percent)
    #print(feh)
    #print(final_prob)
    
    plt.figure(figsize=(7,14))
    plt.subplot(211)
    plt.plot(feh,final_prob,label='Observed Selection metallicity',color = 'blue')
    plt.xlabel('[Fe/H]')
    plt.ylabel('Probability of Observing a Star with our Framework')    
    plt.subplot(212)    
    plt.plot(feh,final_prob_model,label='GALAXIA Metallicity', color = 'red')
    plt.xlabel('[Fe/H]')
    plt.ylabel('Probability of Spawning a Star using GALAXIA')
    plt.legend()
    plt.grid()
    plt.savefig(f'/Users/advik/OPACOS/RESULTS/Framework Inferred {fe_h} and GALAXIA {fe_h} PDF for ({long},{lat}).png')
    #plt.show()

    data_feh = pd.DataFrame(filtered_feh)
    # random number generator

    # Plot pandas histogram from dataframe with df.plot.hist (not df.hist)
    ax = data_feh.plot.hist(bins=30, density=True, edgecolor='red', linewidth=0.125)

    # Save default x-axis limits for final formatting because the pandas kde
    # plot uses much wider limits which usually decreases readability
    #xlim = ax.get_xlim()

    # Plot pandas KDE
    data_feh.plot.density(color='blue', alpha=1, ax=ax) # same as df['var'].plot.kde()

    # Reset x-axis limits and edit legend and add title
    ax.set_xlim([-3.0,0.5])
    ax.legend(labels=['KDE'], frameon=False)
    ax.set_title('Histogram overlaid with KDE', fontsize=14, pad=15)
    plt.savefig(f'/Users/advik/OPACOS/RESULTS/Data {fe_h} Histogram for ({long},{lat}).png')
    #plt.show()
