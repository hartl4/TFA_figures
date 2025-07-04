######################################################################################################
#CODE TO PLOT FIGURE S6
######################################################################################################
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from mpl_toolkits.basemap import shiftgrid
from matplotlib.colors import LogNorm
import cmcrameri.cm as cmc
from mpl_toolkits.axes_grid1 import make_axes_locatable
######################################################################################################
#load datasets
files_tend = [f"/pathto/tend4d_{year}_latlon.nc" for year in range(2020, 2023)] #base tend files
files_surf = [f"/pathto/surf_{year}_latlon.nc" for year in range(2020, 2023)] # base surf files
######################################################################################################
#calculate rainwater concentrations from datasets
for i, (file_tend, file_surf) in enumerate(zip(files_tend, files_surf)):
    # Open the datasets
    ds_tend = xr.open_dataset(file_tend, engine='netcdf4', decode_times=False)
    ds_surf = xr.open_dataset(file_surf, engine='netcdf4', decode_times=False)
    
    # Extract data
    tfa_dep = ds_tend['tfa_wetcnv'].sum(dim='lev') * 1E9  # convert to micrograms
    precip = ds_surf['prec'].sum(dim='lev') * ds_surf['area'].mean(dim='lev') #convert precip from kg/m2 to l
    
    # Convert to numpy arrays
    tfa_dep_np = abs(tfa_dep).values
    precip_np = precip.values
    
    # Exclude values where tfa_dep is 0 or precip is <1 mm/day
    valid_mask = (tfa_dep_np > 0) & (precip_np > 2.4E9)  
    
    # Initialize tfa_conc array with NaNs
    tfa_conc = np.full_like(tfa_dep_np, np.nan)
    
    # Calculate tfa_conc only for valid points
    tfa_conc[valid_mask] = tfa_dep_np[valid_mask] / precip_np[valid_mask]

# Create a masked array instead of extracting values
tfa_conc_masked = np.ma.masked_array(tfa_conc, mask=~valid_mask)

# calc mean over time
tfa_conc_annual = tfa_conc_masked.mean(axis=0) 
vmin, vmax = np.percentile(tfa_conc_annual, [0.5, 99.5])
###############################################################################################
# Create figure and specify axes for map and colorbar
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

# Create the map
m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90, 
            llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax)

# Convert coordinates
lons, lats = np.meshgrid(tfa_dep['lon'], tfa_dep['lat'])
x, y = m(lons, lats)

# Plot
cs = m.contourf(x, y, tfa_conc_annual, cmap=cmc.batlow, 
                norm=LogNorm(vmin=vmin, vmax=vmax))

# Add map features
m.drawcoastlines()
m.drawmapboundary(fill_color='0.3')
m.drawparallels(range(-90, 91, 30), labels=[1, 0, 0, 0], fontsize=7)  # left side only
m.drawmeridians(range(-180, 181, 60), labels=[0, 0, 0, 1], fontsize=7)  # bottom only
# Add labels
#ax.set_title("TFA Annual Mean Precipitation Concentration", fontsize = 10)
ax.set_xlabel('Longitude', fontsize = 10, labelpad=15)
ax.set_ylabel('Latitude', fontsize = 10, labelpad=15)

# Create colorbar 
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
colorbar = plt.colorbar(cs, cax=cax)
colorbar.set_label('Concentration (µg/l)', rotation=270, labelpad=15, fontsize = 10)

plt.tight_layout()
plt.savefig('FigureS6', dpi=300, bbox_inches='tight')
plt.show()