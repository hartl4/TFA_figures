################################################################################################
#CODE TO PLOT FIGURE S8
################################################################################################
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import itertools
import numpy as np

################################################################################################
#LOAD DATASETS
files_tend = [f"/pathto/tend4d_{year}_latlon.nc" for year in range(2000, 2023)]
files_surf = [f"/pathto/surf_{year}_latlon.nc" for year in range(2000, 2023)]

################################################################################################
#CALCULATE ANNUAL AVERAGE RAINWATER CONCENTRATIONS FOR EACH GRIDBOX
tfa_conc_2010 = []
tfa_conc_2022 = []
tfa_conc_11 = []

# Process each year's data
for i, (file_tend, file_surf) in enumerate(zip(files_tend, files_surf)):
    # Open the datasets
    ds_tend = xr.open_dataset(file_tend, engine='netcdf4', decode_times=False)
    ds_surf = xr.open_dataset(file_surf, engine='netcdf4', decode_times=False)
    year = int(file_tend.split('4d_')[1].split('_lat')[0])
    # Extract data
    tfa_dep = ds_tend['tfa_wetcnv'].sum(dim=('time','lev')) * 1E12   #ng
    precip = ds_surf['prec'].sum(dim=('time','lev'))  * ds_surf['area'].mean(dim=('time','lev'))
    
    # Convert to numpy arrays
    tfa_dep_np = abs(tfa_dep).values
    precip_np = precip.values
    
    # Create mask for valid data points
    valid_mask = (tfa_dep_np > 0) & (precip_np > 0) 
    
    # Initialize tfa_conc array with NaNs
    tfa_conc = np.full_like(tfa_dep_np, np.nan)
    
    # Calculate tfa_conc only for valid points
    tfa_conc[valid_mask] = tfa_dep_np[valid_mask] / precip_np[valid_mask]
    tfa_conc_annual =  tfa_conc[valid_mask]

    if year <= 2010:
        tfa_conc_2010.append(tfa_conc_annual)
    elif year == 2011:
        tfa_conc_11.append(tfa_conc_annual)
    else:
        tfa_conc_2022.append(tfa_conc_annual)

tfa_2010_combined = list(itertools.chain.from_iterable(tfa_conc_2010))
tfa_2022_combined = list(itertools.chain.from_iterable(tfa_conc_2022))
tfa_2010_combined_log = np.log10(tfa_2010_combined)
tfa_2022_combined_log = np.log10(tfa_2022_combined )
print(f"Min and max of tfa_2010_combined: {min(tfa_2010_combined_log)}, {max(tfa_2010_combined_log)}")
print(f"Min and max of tfa_2022_combined: {min(tfa_2022_combined_log)}, {max(tfa_2022_combined_log)}")
hist_min = 0.1
hist_max = 200
#################################################################################################################
#plot histogram
# Creating subplots with multiple histograms
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 3), sharey = True)

axes[0].hist(tfa_2010_combined, bins=1000, color='Yellow', edgecolor='black',range=(hist_min, hist_max))
axes[0].set_title('(a) 2000-2010', fontsize = 9)
axes[0].tick_params(axis='both', labelsize=7) 
 
axes[1].hist(tfa_2022_combined, bins=1000, color='Pink', edgecolor='black',range=(hist_min, hist_max))
axes[1].set_title('(b) 2012-2022',fontsize = 9)
axes[1].tick_params(axis='both', labelsize=7) 
 
# Adding labels and title
for ax in axes:
    ax.set_xlabel('TFA concentration (\u03bcg/l)', fontsize = 9)
axes[0].set_ylabel('Frequency',fontsize = 9)
 
plt.tight_layout()
plt.savefig('FigureS8', bbox_inches = 'tight')