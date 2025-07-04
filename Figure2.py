# Description: This code calculates the TFA deposition from variables 'tfa_wetcnv' (wet depostion) and 
#               'tfa_drydep' (dry deposition) from the _tend4d outputs. 
# Author: Lucy Hart
# Date: 09/06/25

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
import seaborn as sns
#IMPORT FILES
#basefiles = f"/INSERTFILEPATH" for year in range(2000, 2023)
#minfiles  = f"/INSERTFILEPATH" for year in range(2000, 2023) 
#maxfiles =[f"/INSERTFILEPATH" for year in range(2000, 2023)] 
#anaestfiles = [f"/INSERTFILEPATH" for year in range(2000, 2023)]
#simplefiles = [f"/INSERTFILEPATH" for year in range(2000, 2023)]
#base42files = [f"/INSERTFILEPATH" for year in range(2000, 2023)]
#surffiles = [f"/INSERTFILEPATH" for year in range(2000, 2023)]

basefiles = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/base/tend4d_{year}_latlon.nc" for year in range(2000, 2023)]# + [f"/home/hartl4/hartl4/outputs/multis_sep/simple99/tend4d_{year}_latlon.nc" for year in range(2002, 2008)] + [f"/home/hartl4/hartl4/outputs/multis_sep/t42simple_07/tend4d_{year}_latlon.nc" for year in range(2008, 2018)]
minfiles  = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/min/tend4d_{year}_latlon.nc" for year in range(2000, 2022)]# + [f"/home/hartl4/hartl4/outputs/multis_sep/simple99/tend4d_{year}_latlon.nc" for year in range(2002, 2008)] + [f"/home/hartl4/hartl4/outputs/multis_sep/t42simple_07/tend4d_{year}_latlon.nc" for year in range(2008, 2018)]
maxfiles = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/max/tend4d_{year}_latlon.nc" for year in range(2000, 2022)]  #+[f"/home/hartl4/hartl4/outputs/multis_sep/simple99/tend4d_{year}_latlon.nc" for year in range(2002, 2008)] + [f"/home/hartl4/hartl4/outputs/multis_sep/t42simple_07/tend4d_{year}_latlon.nc" for year in range(2008, 2018)]
anaestfiles =  (
    [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/anaest/tend4d_{year}_latlon.nc" for year in range(2000, 2023)] 
)
simplefiles = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/simple/tend4d_{year}_latlon.nc" for year in range(2000, 2023)] 
base42files = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/BASE-LOWRES/tend4d_{year}_latlon.nc" for year in range(2000, 2023)] 
surffiles = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/expts/base/surf_{year}_latlon.nc" for year in range(2000, 2023)]

# Initialize a dictionary to hold cumulative and non-cumulative TFA deposition values
tfa_conc = 0
files_list = {'base': basefiles,'min': minfiles, 'max': maxfiles, 'simple': simplefiles, 'anaest': anaestfiles, 'base-lowres': base42files}

# Initialize a dictionary to hold cumulative and non-cumulative TFA deposition values
tfa_cumulative = {'base': {}, 'min' : {}, 'max' : {}, 'simple' : {}, 'anaest' : {}, 'base-lowres' : {}}
tfa_yearly = {'base' : {}, 'min' : {}, 'max': {}, 'simple' : {}, 'anaest' : {}, 'base-lowres' : {}}
#############################################################################################################
#PREP DATA FOR Panels (a) and (b)
#############################################################################################################
#Calculates total deposition (wet + dry) of TFA
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
        year = file.split('4d_')[1].split('_lat')[0]  
        year_int = int(year) 
        prev_year = str(year_int - 1)    
    
    # Add the current year's TFA total to the cumulative deposition
        tfa_cumulative[file_type][year] = tfa_total + tfa_cumulative[file_type].get(prev_year, 0)
    
    # Store the non-cumulative deposition for this year
        tfa_yearly[file_type][year] = tfa_total
#########################################################################################################
#PREP DATA FOR Panels (c) and (d)
#########################################################################################################
tfa_dep_tot = 0
year_dep = 0

# Iterate over the files and calculate cumulative and non-cumulative deposition
for file, surf in zip(basefiles, surffiles):
    # Open the dataset 
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
################################################################################################
#CODE TO PLOT FIGURE
################################################################################################
#generate colorblind friendly colormap

cmap = cmc.batlow
n_sources = len(files_list)
colors = sns.color_palette("colorblind", n_sources)


fig, ax = plt.subplots(2,2, sharey=False, figsize=(10, 5), gridspec_kw={'width_ratios': [1, 1.5]})

# Plot LHS: Time Series Graphs
for i, file_type in enumerate(files_list.keys()):
    years = sorted(tfa_cumulative[file_type].keys())
    mean_cumulative_values = [tfa_cumulative[file_type][year] for year in years]
    mean_yearly_values = [tfa_yearly[file_type][year] for year in years]

    # Plot cumulative and non-cumulative deposition with unique color
    ax[1, 0].plot(years, mean_cumulative_values, marker='',linestyle='--',  label=f'{file_type.upper()}', color=colors[i])
    ax[0, 0].plot(years, mean_yearly_values, marker='', linestyle='--', label=f'{file_type.upper()}', color=colors[i])

# Collect all years from the data
all_years = sorted(set().union(*[tfa_cumulative[ft].keys() for ft in files_list.keys()]))

# Add labels and legend
lon = ds['lon']
lat = ds['lat']
ax[1,0].set_title('(b) Cumulative Deposition', fontsize = 10)
ax[0,0].set_title('(a) Annual Deposition', fontsize = 10)
ax[1,0].set_xlabel('Year', fontsize = 10)
ax[0,0].set_ylabel('TFA Deposition (Gg/yr)', fontsize = 10)
ax[1,0].set_ylabel('TFA Deposition (Gg)', fontsize = 10)
ax[0,0].legend(loc='upper left', prop={'size': 7})
ax[0,0].grid(True)
ax[1,0].grid(True)
ax[0,0].set_xticks(years[::2])
ax[1,0].set_xticks(years[::2])
ax[0,0].set_xticklabels([]) 
ax[1,0].set_xticklabels(years[::2], rotation=45)

#PLOT DEPOSITION MAP
lon, lat = np.meshgrid(tfa_dep['lon'], tfa_dep['lat'])
contour = ax[0, 1].contourf(lon, lat, tfa_dep.values, cmap=cmc.batlow)
cbar = fig.colorbar(contour, ax=ax[0, 1], orientation='vertical')
ax[0,1].set_title("(c) Average TFA Deposition 2000-2022 (µg/m\u00b2/yr)",fontsize = 10 )
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax[0, 1])
m.drawmapboundary(fill_color='0.3')
m.drawcoastlines()

#PLOT CUMULATIVE DEPOSITION MAP
lon2, lat2 = np.meshgrid(tfa_dep_tot['lon'], tfa_dep_tot['lat'])
contour2 = ax[1, 1].contourf(lon2, lat2, tfa_dep_tot.values, cmap=cmc.batlow)
cbar = fig.colorbar(contour2, ax=ax[1, 1], orientation='vertical')
ax[1,1].set_title("(d) Cumulative TFA Deposition 2000-2022 (µg/m\u00b2)",fontsize = 10 )
m2 = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax[1, 1])
m2.drawmapboundary(fill_color='0.3')
m2.drawcoastlines()

# Adjust layout
plt.tight_layout()
plt.savefig('Figure2', bbox_inches = 'tight')
plt.show()
