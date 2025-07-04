################################################################################################
#CODE TO PLOT FIGURE S4
################################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cmcrameri.cm as cmc
###################################################################################################
#CALCULATE PRODUCTION IN EACH REGION
# Constants
MW_air = 28.97
MW_TFAC = 132.4
MW_TFF = 116
MW_CF3CHO = 98
MW_TFA = 114

# Initialize dictionaries for each region
regions_prod = {
    'North Pole (>66°N)': {'lat_range': (66, 91), 'monthly_data': {}},
    'Midlat NH (23.27-66°N)': {'lat_range': (23.27, 66), 'monthly_data': {}},
    'Tropics NH (0-23.27°N)': {'lat_range': (0, 23.27), 'monthly_data': {}},
    'Tropics SH (23.27°S-0)': {'lat_range': (-23.27, 0), 'monthly_data': {}},
    'Midlat SH (66-23.27°S)': {'lat_range': (-66, -23.27), 'monthly_data': {}},
    'South Pole (<66°S)': {'lat_range': (-90, -66), 'monthly_data': {}}
}
# Initialize monthly data for each region
for region_p in regions_prod.values():
    region_p['monthly_data'] = {str(month).zfill(2): [] for month in range(1, 13)}

# File paths
files = [f"/pathto/chem_{year}_latlon.nc" for year in range(2000, 2022)] # base chem files
files_surf = [f"/pathto/surf_{year}_latlon.nc" for year in range(2000, 2022)] # base surf files

# calculate TFA production from chem files
for file, file_surf in zip(files, files_surf):
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    ds2 = xr.open_dataset(file_surf, engine='netcdf4', decode_times=False)
    
    # Process each region
    for region_name, region_prod in regions_prod.items():
        lat_min, lat_max = region_prod['lat_range']
        
        # Select region data
        ds_region = ds.sel(lat=slice(lat_min, lat_max))
        ds2_region = ds2.sel(lat=slice(lat_min, lat_max))
        
        areax = ds2_region['area']
        area_sum = areax.sum()
        
        # Calculate TFA for each month in the region
        for month in range(1, 13):
            month_idx = month - 1
            
            tfa_chem = (
                (np.sum(ds_region['chem193'].isel(time=month_idx)) * MW_TFA/MW_air * 3600) +
                (np.sum(ds_region['chem192'].isel(time=month_idx)) * MW_TFA/MW_air * 3600) +
                (np.sum(ds_region['chem191'].isel(time=month_idx)) * MW_TFA/MW_air * 3600)
            )
            
            tfa_total_prod = abs(np.sum(tfa_chem) * 1e-6) #kg to Gg
            tfa_byarea = tfa_total_prod * 1E15/area_sum # Gg to micrograms/m2
            
            region_prod['monthly_data'][str(month).zfill(2)].append(tfa_byarea.item())

# Calculate monthly averages for each region
for region_prod in regions_prod.values():
    region_prod['monthly_averages'] = {
        month: np.mean(values) 
        for month, values in region_prod['monthly_data'].items()
    }
##################################################################################################
#CALCULATE DEPOSITION IN EACH REGION
# Initialize dictionaries for each region
regions = {
    'North Pole (>66°N)': {'lat_range': (66, 91), 'monthly_data': {}},
    'Midlat NH (23.27-66°N)': {'lat_range': (23.27, 66), 'monthly_data': {}},
    'Tropics NH (0-23.27°N)': {'lat_range': (0, 23.27), 'monthly_data': {}},
    'Tropics SH (23.27°S-0)': {'lat_range': (-23.27, 0), 'monthly_data': {}},
    'Midlat SH (66-23.27°S)': {'lat_range': (-66, -23.27), 'monthly_data': {}},
    'South Pole (<66°S)': {'lat_range': (-90, -66), 'monthly_data': {}}
}

# Initialize monthly data for each region
for region in regions.values():
    region['monthly_data'] = {str(month).zfill(2): [] for month in range(1, 13)}

# File paths
files = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/base/tend4d_{year}_latlon.nc" for year in range(2000, 2022)]
files_surf = [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/expts/base/surf_{year}_latlon.nc" for year in range(2000, 2022)]

# Process each year's data
for file, file_surf in zip(files, files_surf):
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    ds2 = xr.open_dataset(file_surf, engine='netcdf4', decode_times=False)
    
    # Process each region
    for region_name, region_dep in regions.items():
        lat_min, lat_max = region_dep['lat_range']
        
        # Select region data
        ds_region = ds.sel(lat=slice(lat_min, lat_max))
        ds2_region = ds2.sel(lat=slice(lat_min, lat_max))
        
        areax = ds2_region['area']
        area_sum = areax.sum()
        
        # Calculate TFA for each month in the region
        for month in range(1, 13):
            month_idx = month - 1
            
            tfa_dep = (
                (np.sum(ds_region['tfa_wetcnv'].isel(time=month_idx)))  +
                (np.sum(ds_region['tfa_drydep'].isel(time=month_idx)))  
            )
            
            tfa_total_dep = abs((tfa_dep) * 1e-6)
            tfa_byarea = np.divide(tfa_total_dep * 1E15,area_sum)
            
            region_dep['monthly_data'][str(month).zfill(2)].append(tfa_byarea.item())

# Calculate monthly averages for each region
for region_dep in regions.values():
    region_dep['monthly_averages'] = {
        month: np.mean(values) 
        for month, values in region_dep['monthly_data'].items()
    }

##################################################################################################
#PLOTTING
#PLOT SUBPLOTS
cmap = cmc.batlow
n_sources = len(regions_prod)
colors = [cmap(i / n_sources) for i in range(n_sources)]

fig, ax = plt.subplots(2, sharey=True, sharex=True, figsize=(10, 8))

#plot prod data on ax1
for (region_name, region_prod), color in zip(regions_prod.items(), colors):
    months_prod = list(region_prod['monthly_averages'].keys())
    averages_prod = list(region_prod['monthly_averages'].values())
    
    # Plot line and points
    ax[0].plot(months_prod, averages_prod, marker='o', linestyle='-', linewidth=2, 
             label=region_name, color=color)
    
    # Calculate and plot standard deviation
    stds = [np.std(region_prod['monthly_data'][month]) for month in months_prod]


#plot dep data on ax2
for (region_name, region_dep), color in zip(regions.items(), colors):
    months_dep = list(region_dep['monthly_averages'].keys())
    averages_dep = list(region_dep['monthly_averages'].values())
    
    # Plot line and points
    ax[1].plot(months_dep, averages_dep, marker='o', linestyle='-', linewidth=2, 
        label=region_name, color=color)
    
    # Calculate and plot standard deviation
    stds = [np.std(region_dep['monthly_data'][month]) for month in months_dep]


ax[1].set_xlabel('Month')
ax[0].set_title('(a) TFA Production', fontsize = 10)
ax[1].set_title('(b) TFA Deposition', fontsize = 10)
ax[0].set_ylabel('Average TFA Production (µg/m\u00b2)', fontsize = 10)
ax[1].set_ylabel('Average TFA Depositionion (µg/m\u00b2)',fontsize = 10)

# Customize x-axis
month_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
ax[1].set_xticks(months_prod, month_labels)

# Add legend
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=7, frameon=False,bbox_to_anchor=(0.5, -0.00))

plt.savefig('FigureS4.png')
plt.show()