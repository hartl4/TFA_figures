#IMPORT MODULES 
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmcrameri.cm as cmc

#IMPORT FILES
basefiles = f"/INSERTFILEPATH" for year in range(2000, 2023)
minfiles  = f"/INSERTFILEPATH" for year in range(2000, 2023) 
maxfiles =[f"/INSERTFILEPATH" for year in range(2000, 2023)] 
anaestfiles = [f"/INSERTFILEPATH" for year in range(2000, 2023)]
simplefiles = [f"/INSERTFILEPATH" for year in range(2000, 2023)]
base42files = [f"/INSERTFILEPATH" for year in range(2000, 2023)]
surffiles = [f"/INSERTFILEPATH" for year in range(2000, 2023)] #_surf files

# Initialize a dictionary to hold cumulative and non-cumulative TFA deposition values
tfa_conc = 0
files_list = {'base': basefiles,'min': minfiles, 'max': maxfiles, 'simple': simplefiles, 'anaest': anaestfiles, 'base-lowres': base42files}

# Initialize a dictionary to hold cumulative and non-cumulative TFA deposition values
tfa_cumulative = {'base': {}, 'min' : {}, 'max' : {}, 'simple' : {}, 'anaest' : {}, 'base-lowres' : {}}
tfa_yearly = {'base' : {}, 'min' : {}, 'max': {}, 'simple' : {}, 'anaest' : {}, 'base-lowres' : {}}
#########################################################################################################################################################
#PREP DATA FOR Panels (a) and (b)
# Iterate over the files and calculate cumulative and non-cumulative deposition
for file_type, files in files_list.items():
    for file in files:
    # Open the dataset using xarray and decode times as False to avoid time conversion
        ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
        year = file.split('4d_')[1].split('_lat')[0]      
        year_int = int(year)  
    
    # Extract the variables for dry and wet TFA deposition
        tfa_dry = ds['tfa_drydep'] #kg/gridbox/month
        tfa_wet = ds['tfa_wetcnv'] #kg/gridbox/month

    # Calculate total TFA deposition (dry + wet) for the current file (yearly total)
        tfa_total = abs(np.sum(tfa_dry) * 1e-6) + abs(np.sum(tfa_wet) * 1e-6)  # Convert to Gg
        year = file.split('4d_')[1].split('_lat')[0]      # Extract the year from the filename (assuming consistent naming convention)
        year_int = int(year)  # Convert year to integer for manipulation
        prev_year = str(year_int - 1)     # Get the previous year as a string
    
    # Add the current year's TFA total to the cumulative deposition
        tfa_cumulative[file_type][year] = tfa_total + tfa_cumulative[file_type].get(prev_year, 0)
    
    # Store the non-cumulative deposition for this year
        tfa_yearly[file_type][year] = tfa_total
#########################################################################################################################################################
#PREP DATA FOR Panels (c) and (d)
tfa_dep_tot = 0
year_dep = 0

# Iterate over the files and calculate cumulative and non-cumulative deposition
for file, surf in zip(basefiles, surffiles):
    # Open the dataset using xarray and decode times as False to avoid time conversion
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    sa = xr.open_dataset(surf, engine='netcdf4', decode_times=False)
    # Extract the variables for dry and wet TFA deposition
    tfa_dry = ds['tfa_drydep']#kg/gridbox/month
    tfa_wet = ds['tfa_wetcnv']#kg/gridbox/month
    area = sa['area'].mean(dim= ('time','lev')) #gridbox area in m2
    tfa_total = abs(tfa_dry.sum(dim=('time','lev')) + tfa_wet.sum(dim=('time','lev'))) #calculate annual tfa deposition for each lat/lon
    tfa_total_byarea = np.divide(tfa_total, area)
    tfa_dep_tot += tfa_total_byarea*1E9
    year_dep += 1

tfa_dep = tfa_dep_tot/year_dep  #average annual dep
#########################################################################################################################################################
#CODE TO PLOT FIGURE
fig, ax = plt.subplots(2,2, sharey=False, figsize=(15, 7), gridspec_kw={'width_ratios': [1, 1.5]})

#PLOT LHS AS THE LINE GRAPHS
for file_type in files_list.keys():
    years = sorted(tfa_cumulative[file_type].keys())
    mean_cumulative_values = [tfa_cumulative[file_type][year] for year in years]
    mean_yearly_values = [tfa_yearly[file_type][year] for year in years]
    # Choose marker style and size based on file_type
# Choose marker style based on file_type
    if file_type == 'basefiles':
        marker_style = 'o'
    elif file_type == 'base42files':
        marker_style = 's'
    elif file_type == 'anaestfiles':
        marker_style = '^'
    else:
        marker_style = 'o'
   
    # Plot cumulative and non-cumulative deposition
    ax[1, 0].plot(years, mean_cumulative_values, marker=marker_style, label=f'{file_type.upper()}')
    ax[0, 0].plot(years, mean_yearly_values, marker=marker_style, linestyle='--', label=f'{file_type.upper()}')

# Add labels and legend
lon = ds['lon']
lat = ds['lat']
ax[1,0].set_title('(b) Cumulative Deposition')
ax[0,0].set_title('(a) Annual Deposition')
ax[1,0].set_xlabel('Year')
ax[0,0].set_ylabel('TFA Deposition (Gg/yr)')
ax[1,0].set_ylabel('TFA Deposition (Gg)')
ax[0,0].legend(loc='upper left')
ax[0,0].grid()
ax[1,0].grid()
#ax[0,0].tick_params(axis='x', rotation=45)
ax[0,0].set_xticklabels([]) 
ax[1,0].tick_params(axis='x', rotation=45) 

#PLOT DEPOSITION MAP
lon, lat = np.meshgrid(tfa_dep['lon'], tfa_dep['lat'])
contour = ax[0, 1].contourf(lon, lat, tfa_dep.values, cmap=cmc.batlow)

#cbar = fig.colorbar(contour, ax=ax[0, 1], orientation='vertical', pad=0.01)
cbar = fig.colorbar(contour, ax=ax[0, 1], orientation='vertical')
ax[0,1].set_title("(c) Average TFA Deposition 2000-2022 (µg/m\u00b2/yr)" )
 # Adjust tick density for clarity
#ax[0,1].set_title('(d) TFA Deposition (kg/yr)', rotation=0, labelpad=0)
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax[0, 1])
m.drawmapboundary(fill_color='0.3')
m.drawcoastlines()

#PLOT CUMULATIVE DEPOSITION
lon2, lat2 = np.meshgrid(tfa_dep_tot['lon'], tfa_dep_tot['lat'])
contour2 = ax[1, 1].contourf(lon2, lat2, tfa_dep_tot.values, cmap=cmc.batlow)

cbar = fig.colorbar(contour2, ax=ax[1, 1], orientation='vertical')
ax[1,1].set_title("(d) Cumulative TFA Deposition 2000-2022 (µg/m\u00b2)" )
 # Adjust tick density for clarity
#ax[0,1].set_title('(d) TFA Deposition (kg/yr)', rotation=0, labelpad=0)
m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax[1, 1])
m2.drawmapboundary(fill_color='0.3')
m2.drawcoastlines()

# Adjust layout
plt.tight_layout()
plt.savefig('grl2', bbox_inches = 'tight')
plt.show()

