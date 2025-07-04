# Description: This code uses chemical fluxes (kg/gridbox/month) 
#               from the _chem model outputs to calculate tge TFA yield for each source gas,
#               as described in Text S3
# Author: Lucy Hart
# Date: 09/06/25

#IMPORT MODULES
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
from netCDF4 import Dataset
import cmcrameri.cm as cmc

#files
#files_base = [f"/INSERTFILEPATH" for year in range(2000, 2023)] # tend4d files for BASE
#files_chembase = [f"/INSERTFILEPATH" for year in range(2000, 2023)] #CHEM files for BASE
#files_cheman = [f"INSERTFILEPATH" for year in range(2000, 2023)] #CHEM files for BASE-AN
#files_anaest = [f"/INSERTFILEPATH"] for year in range(2000, 2023)] #CHEM files for BASE-AN 
#files_simple = [f"/INSERTFILEPATH"] for year in range(2000, 2023)] #CHEM files for SIMPLE-Chem 
files_base = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/base/tend4d_{year}_latlon.nc" for year in range(2000, 2023)]
files_chembase = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/base/chem_{year}_latlon.nc" for year in range(2000, 2023)]
files_cheman = (
    [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/anaest/chem_{year}_latlon.nc" for year in range(2000, 2023)]
)
files_simple = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/simple/chem_{year}_latlon.nc" for year in range(2000, 2023)]
files_anaest =  (
    [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/anaest/tend4d_{year}_latlon.nc" for year in range(2000, 2023)]
)

##########################################################################################
#CODE TO CALCULATE PRODUCTION OF TFA FROM EACH SOURCE GAS
##########################################################################################
#load dicts
#EXPTS/BASE
my_dict = {
    'HCFC-123': {"MW": 152.9E-3, "prod": ['chem078', 'chem182'], "loss": 'chem077', "prec": 132.4E-3, "yield":[], "which": 'tfac'},
    'HCFC-124': {"MW": 136.5E-3, "prod": ['chem082', 'chem184'], "loss": 'chem081', "prec": 116E-3, "yield":[], "which": 'tff'},
    'HCFC-133a': {"MW": 118.5E-3, "prod": ['chem086', 'chem186'] , "loss": 'chem085', "prec": 132.4E-3, "yield":[], "which": 'tfac'},
    'HFC-134a': {"MW": 102.0E-3, "prod": ['chem091', 'chem190'], "loss": 'chem090', "prec": 116E-3, "yield":[], "which": 'tff'}}
my_dict2 = {
    'tff': {"MW": 116E-3, "prod": ['chem082', 'chem184', 'chem091', 'chem188', 'chem097', 'chem190', 'chem117' ], "loss": 'chem199',  "yield":[]},
    'tfac': {"MW": 132.4E-3, "prod": ['chem078', 'chem182', 'chem086', 'chem186'], "loss": 'chem198',  "yield":[]},
    'cf3cho': {"MW": 98E-3, "prod": ['chem102', 'chem107',  'chem192', 'chem194', 'chem196'], "loss": 'chem200',  "yield":[]},
    'isof': {"MW": 164E-3, "prod": ['chem126'], "loss": 'chem201',  "yield":[]},
    'tfa': {"MW": 114E-3, "prod": [], "loss": [],  "yield":[]} }

# MWs IN g/mol
MW_AIR = 28.97 
MW_TFA = 114

#initialise variables for yields calc
output_data = []
all_losses = {}
all_losses_source = {}
all_prod = {}
all_prod_source = {}
total_loss_sum = 0
total_loss_sum_source = 0
total_prod_sum = 0
total_prod_sum_source = 0

# Initialize collections to store data across multiple years
years_data = []
tfa_production_by_year = {}

#calculate yields and production
for file, file_tend in zip(files_cheman, files_anaest):
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    year = file.split('em_')[1].split('_lat')[0]
    year_data = {"year": year, "tfa_yield": {}}
    year_data_tfa_prod = {"year": year, "tfa_production_Gg": {}}  # New dictionary for TFA production in Gg
    
    # First pass: Calculate TFA yields for each precursor
    for key2 in my_dict2:
        total_loss_sum = 0
        total_prod_sum = 0
        
        # Convert to lists
        losses = my_dict2[key2]["loss"]
        prods = my_dict2[key2]["prod"]
        
        if not isinstance(losses, list):
            losses = [losses]
        if not isinstance(prods, list):
            prods = [prods]
        
        # Store in dictionaries
        all_losses[key2] = losses
        all_prod[key2] = prods
        
        # Calculate chemical flux in kg/year
        for loss in losses:
            total_loss_sum += (ds[loss].sum()*3600*MW_TFA/MW_AIR)
        for prod in prods:
            total_prod_sum += (ds[prod].sum()*3600*my_dict2[key2]["MW"]/MW_AIR).values
        
        # Calculate TFA yield for the current key
        if total_loss_sum != 0:
            tfa_yield = (total_loss_sum /MW_TFA) /(total_prod_sum / my_dict2[key2]["MW"])
            year_data["tfa_yield"][key2] = tfa_yield
        else:
            year_data["tfa_yield"][key2] = 0.0
    
    # Second pass: Calculate TFA production in Gg from each source gas
    for key in my_dict:
        total_prod_sum_source = 0
        
        # Get the precursor key from my_dict for this source gas
        precursor_key = my_dict[key]["which"]
        prods_source = my_dict[key]["prod"]
        
        if not isinstance(prods_source, list):
            prods_source = [prods_source]
        
        # Store in dictionaries
        all_prod_source[key] = prods_source
        
        # Calculate production in kg/year
        for prod_source in prods_source:
            total_prod_sum_source += (ds[prod_source].sum()*3600*my_dict[key]["prec"]/MW_AIR).values
        
        # Get the TFA yield for this precursor
        if precursor_key in year_data["tfa_yield"]:
            precursor_tfa_yield = year_data["tfa_yield"][precursor_key]
            
            # Calculate TFA production in Gg
            tfa_production = precursor_tfa_yield * total_prod_sum_source
            tfa_production_gg= tfa_production*1E-6*(MW_TFA/my_dict[key]["prec"])

            # Store the year
            years_data.append(year)
    
            # Store TFA production data for this year
            if year not in tfa_production_by_year:
                tfa_production_by_year[year] = {}

            tfa_production_by_year[year][key] = tfa_production_gg  # Store each source gas separately
            
        else:
            year_data_tfa_prod["tfa_production_Gg"][key] = 0.0


# tfa production line
#calculates production of tfa from precursor hydrolysis
years = sorted(tfa_production_by_year.keys())
for file, file_tend in zip(files_cheman, files_anaest):
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    ds_tend = xr.open_dataset(file_tend, engine='netcdf4', decode_times=False)
    #tfa_prod_terms = np.sum(ds['chem191'] )+ np.sum(ds['chem192']) + np.sum(ds['chem193']) + np.sum(ds['chem119'])  #BASE run
    tfa_prod_terms = np.sum(ds['chem200']) + np.sum(ds['chem199']) + np.sum(ds['chem201']) + np.sum(ds['chem198']) + np.sum(ds['chem119']) + np.sum(ds['chem129']) #BASE-AN run
    tfa_prod_total = (tfa_prod_terms)*3600 *MW_TFA/MW_AIR *1E-6
    year = file.split('em_')[1].split('_lat')[0]
    print(tfa_prod_total.shape)
    tfa_production[year] = tfa_prod_total.item()
tfa = [tfa_production[year] for year in years]

#prepare production data for plotting
sources = set()

data = {source: [] for year in years for source in tfa_production_by_year[year].keys()}

# Collect data for each source
for year in years:
    for source, value in tfa_production_by_year[year].items():
        data[source].append(value)
    sources.update(tfa_production_by_year[year].keys())

sources = sorted(sources)  # Ensure consistent ordering

# Convert to numpy array for stacked plotting
y_values = np.array([data[source] for source in sources])
#############################################################################################
#PANEL2 DATA PREP
#############################################################################################
tfa_base = 0 
tfa_simple = 0
tfa_anaest = 0
years_2 = 0
years_a = 0

# calculate production sum across chemicals
for file_chembase in files_chembase:
    year_file = int(file_chembase.split('em_')[1].split('_lat')[0])
    ds =  xr.open_dataset(file_chembase, engine='netcdf4', decode_times=False)
    prod_base = ds['chem191'] + ds['chem192'] + ds['chem193']
    tfa_base += prod_base.sum(axis=(0,1,3)) * 3600 * MW_TFA/MW_AIR * 1E-6  # Gg units
    years_2 += 1
    
for file_simple in files_simple:
    ds = xr.open_dataset(file_simple, engine='netcdf4', decode_times=False)
    prod_simple = ds['chem130'] + ds['chem129'] + ds['chem128']
    tfa_simple += prod_simple.sum(axis=(0,1,3)) * 3600 * MW_TFA/MW_AIR * 1E-6  # Gg units

for file_anaest in files_cheman:
    ds = xr.open_dataset(file_anaest, engine='netcdf4', decode_times=False)
    prod_anaest =ds['chem200'] + ds['chem199'] + ds['chem201'] + ds['chem198'] + ds['chem119'] + ds['chem129']
    tfa_anaest += prod_anaest.sum(axis=(0,1,3)) * 3600 * MW_TFA/MW_AIR * 1E-6  # Gg units
    years_a += 1

#annual mean across time period for BASE, SIMPLE, and BASE-AN runs
tfa_base_mean = tfa_base/years_a
tfa_simple_mean = tfa_simple/years_a
tfa_anaest_mean = tfa_anaest/years_a

latitudes = ds['lat']

#############################################################################################
#PLOTTING
#############################################################################################
#generate colorblind friendly colormap
cmap = plt.get_cmap('viridis')
n_sources = len(sources)+1
colors = [cmap(i / n_sources) for i in range(n_sources)]

# Plot the stacked area chart
fig, ax = plt.subplots(2, sharey=False, figsize=(8, 9))

# Plotting individual components
ax[0].plot(years, tfa, label='Total TFA ', color='black', zorder = 3) #plot total tfa production line
stack = ax[0].stackplot(years, y_values, labels=sources, alpha=0.7, zorder=1, colors=colors)
stack_top = np.sum(y_values, axis=0) 
#Fill between stack top and black line
ax[0].fill_between(years, stack_top, tfa, where=(tfa > stack_top), 
                   interpolate=True, color=colors, alpha=0.4, label='Other', zorder=2)
# Labels and title
ax[0].set_xticks(years[::2])
ax[0].set_ylabel("TFA Production (Gg/yr)", fontsize=10)
ax[0].set_xlabel('Year', fontsize=10)
ax[0].set_title("(a) TFA Production by Source (BASE-AN)",fontsize=10)
ax[0].legend(loc="upper left", bbox_to_anchor=(0, 1),  prop={'size': 7})

#AX2 PLOT
ax[1].plot(latitudes, tfa_base_mean, label = 'BASE', linestyle='-')
ax[1].plot(latitudes, tfa_simple_mean, label = 'SIMPLE-chem', linestyle='-')
ax[1].plot(latitudes, tfa_anaest_mean, label = 'BASE-AN', linestyle='-')
ax[1].set_xticks(range(-90, 91, 15))
ax[1].set_title("(b) Zonal and Annual Mean TFA Production (2000-2022)", fontsize=10)
ax[1].set_xlabel('Latitude (°)',fontsize=10)
ax[1].set_ylabel('TFA Production (Gg/yr)',fontsize=10)
ax[1].legend(loc="upper left", bbox_to_anchor=(0, 1) ,prop={'size': 7})

plt.savefig('Figure1', bbox_inches='tight')
plt.show()
