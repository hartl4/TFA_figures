######################################################################################################
#CODE TO PLOT FIGURE S5
######################################################################################################
import matplotlib.path as mpath
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
from mpl_toolkits.basemap import shiftgrid
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs, cartopy.feature as cfeature
from scipy.interpolate import griddata
import cmcrameri.cm as cmc

######################################################################################################
#import files 
files_tend63 = [f"/pathto/tend4d_{year}_latlon.nc" for year in range(2000, 2023)] #BASE
surf63 = "/pathto/base/surf_2010_latlon.nc" 
surf63= f"/pathto/surf_2010_latlon.nc" 
sa63 =xr.open_dataset(surf63,engine='netcdf4',decode_times=False)
files_tend42 = [f"/pathto/base42/tend4d_{year}_latlon.nc" for year in range(2000, 2023)] #BASE LOW-RES
#surf42 = f"/pathto/base42/surf_2010_latlon.nc" 
surf42= f"/pathto/surf_2010_latlon.nc" 
sa42 =xr.open_dataset(surf42,engine='netcdf4',decode_times=False)
######################################################################################################

# Process TFA deposition data
def process_and_plot(files, ax, area_data, title, add_colorbar=True):

    tfa_dep_tot = None  
    
    for file in files:
        ds = xr.open_dataset(file)
        # Standardize coordinates
        if 'lat' in ds.coords:
            ds = ds.assign_coords(lat=ds.lat.astype('float64').round(4))
        if 'lon' in ds.coords:
            ds = ds.assign_coords(lon=ds.lon.astype('float64').round(4))
        
        data_to_add = abs(ds['tfa_wetcnv'] + ds['tfa_drydep']).sum(dim=('time', 'lev')) 
        
        if tfa_dep_tot is None:
            tfa_dep_tot = data_to_add
        else:
            # Subsequent iterations: align and add
            tfa_dep_tot, data_to_add = xr.align(tfa_dep_tot, data_to_add, join='override')
            tfa_dep_tot += data_to_add
    tfa_dep = tfa_dep_tot/len(files)
        
        # Get the area data
    if isinstance(area_data, xr.DataArray):
        area = area_data[0, 0, :, :]
            
        # Calculate TFA - handle potential grid mismatches
    if tfa_dep.shape == area.shape:
        tfa = (tfa_dep * 1E9) / area
        
        # Get the latitude data
    lat_data = ds['lat']
        
        # Apply latitude filter and plot
    if hasattr(lat_data, 'ndim') and lat_data.ndim == 2:
            # For 2D lat grids
        mask = lat_data > lat_lims[0]
        tfa_masked = tfa.where(mask)
        
    else:
            # For 1D lat coordinate
        if 'lat' in tfa.dims:
                # If latitude is a dimension, use where with a condition on that dimension
            mask = lat_data > lat_lims[0]
            tfa_masked = tfa.where(lambda x: x.lat > lat_lims[0])

        else:
            print("Converting to 2D grid for plotting")
                # Extract coordinates
            lats = lat_data.values if hasattr(lat_data, 'values') else lat_data
            lons = ds['lon'].values if hasattr(ds['lon'], 'values') else ds['lon']
            values = tfa.values
        
            if lats.ndim == 1 and lons.ndim == 1 and values.ndim == 2:
                    # Create 2D meshgrid
                lon_grid, lat_grid = np.meshgrid(lons, lats)
                    
                    # Create mask
                mask_indices = np.where(lats > lat_lims[0])[0]
                masked_values = values.copy()
                    
                    # Set values outside the desired latitude range to NaN
                for i in range(len(lats)):
                    if i not in mask_indices:
                        masked_values[i, :] = np.nan
                    
                    # Plot using pcolormesh
                im = ax.pcolormesh(lon_grid, lat_grid, masked_values, 
                                      transform=ccrs.PlateCarree())
                return im

        # Use xarray's built-in plotting with control over colorbar
        im = tfa_masked.plot(ax=ax, transform=ccrs.PlateCarree(), add_colorbar=add_colorbar, cmap = cmc.batlow)
        return im
        


######################################################################################################
# Create a figure with two subplots
fig = plt.figure(figsize=(16, 8))
# Use gridspec for better control of layout with shared colorbar
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.05)
ax1 = fig.add_subplot(gs[0, 0], projection=ccrs.NorthPolarStereo())
ax2 = fig.add_subplot(gs[0, 1], projection=ccrs.NorthPolarStereo())

lat_lims = [70, 90]
#####################################################################################
#Function to compute circle boundary for map
def polarCentral_set_latlim(lat_lims, ax):
    ax.set_extent([-180, 180, lat_lims[0], lat_lims[1]], ccrs.PlateCarree())
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
###################################################################################
#PLOTTING
def sp_map(*nrs, projection = ccrs.PlateCarree(), **kwargs):
    return plt.subplots(*nrs, subplot_kw={'projection':projection}, **kwargs)
def add_map_features(ax):
    ax.coastlines()
    gl = ax.gridlines(color='black')
    ax.add_feature(cfeature.BORDERS);
    gl = ax.gridlines()
    gl.xlabels_top = False
    gl.ylabels_right = False

# Process and plot for ds1 (left plot)
im1 = process_and_plot(files_tend63, ax1, sa63['area'], "BASE run", add_colorbar=False)
polarCentral_set_latlim(lat_lims, ax1)
add_map_features(ax1)
ax1.set_title('(a) BASE run', fontsize=14)

# Process and plot for ds2 (right plot)
im2 = process_and_plot(files_tend42, ax2, sa42['area'], "BASE-LOWRES run", add_colorbar=False)
polarCentral_set_latlim(lat_lims, ax2)
add_map_features(ax2)
ax2.set_title('(b) BASE-LOWRES run', fontsize=14)

# ADD ANNOTATIONS FOR SITE LOCATIONS 
site_locations = {
    'Mt Oxford Icefield': {'lat': 82.19, 'lon': -72.96},
    'Devon Icecap': {'lat': 75.2, 'lon': -82.7},
    'Lomonosovfonna': {'lat': 78, 'lon': 17}
}

offsets = {
    'Mt Oxford Icefield': {'lat_offset': 0, 'lon_offset': 10},
    'Devon Icecap': {'lat_offset': -2, 'lon_offset': -10},
    'Lomonosovfonna': {'lat_offset': 2, 'lon_offset': -0.3}
}

for ax in [ax1, ax2]:
    for site_name, coords in site_locations.items():
        lat_offset = offsets[site_name]['lat_offset']
        lon_offset = offsets[site_name]['lon_offset']
        #add text lables
        ax.text(coords['lon'] + lon_offset, coords['lat'] + lat_offset, site_name,
                transform=ccrs.PlateCarree(), fontsize=12, color='black',
                ha='center', va='center', fontweight='bold',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle="round,pad=0.3"))
        # Add dots at exact site locations
        ax.scatter(coords['lon'], coords['lat'], color='white', s=50,
                  transform=ccrs.PlateCarree(), edgecolor='black', zorder=5)

# Add a single colorbar for both plots
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
cbar = fig.colorbar(im1, cax=cbar_ax)
cbar.set_label('TFA Deposition (µg/m\u00b2/yr)', fontsize=12)

# Adjust layout and add a main title
plt.tight_layout(rect=[0, 0, 0.9, 0.95])  # Make room for colorbar and title
plt.savefig('FigureS5', bbox_inches = 'tight')
plt.show()