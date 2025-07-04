##########################################################################
#CODE TO PLOT FIGURE S2
##########################################################################

#IMPORT MODULES
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cftime
import calendar

##########################################################################
species_list = ['hfc134a', 'hfc143a','hfc227ea', 'hfc245fa', 'hfc365mfc', 'hfc4310mee']
files = [f"pathto/{species}_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-0-3-0_gr1z_0001-2022.nc" for species in species_list]

# Create figure for plotting
fig, axs = plt.subplots(3, 2, figsize=(9, 7), sharex=True)
axs = axs.flatten()

# Create a list to store dataframes for each species
all_data = []

for i, (species, file) in enumerate(zip(species_list, files)):
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    da = ds[species]  # Variable has the same name as species
    
    time = xr.coding.times.decode_cf_datetime(
        ds['time'].values,
        units=ds['time'].attrs['units'],
        calendar=ds['time'].attrs.get('calendar', 'standard')
    )
    
    # Convert cftime objects to decimal years
    decimal_years = np.array([
        t.year + (t.timetuple().tm_yday - 1) / (366 if calendar.isleap(t.year) else 365)
        for t in time
    ])
    da = da.assign_coords(time=decimal_years)
    
    # Filter for years from 2000-2022
    year_mask = (decimal_years >= 2000) & (decimal_years < 2023)
    da = da.isel(time=year_mask)
    
    # Calculate hemispheric and global mean
    weights = np.cos(np.deg2rad(da.lat))
    global_mean = da.weighted(weights).mean(dim='lat')
    nh_mean = da.sel(lat=slice(0, 90)).weighted(weights.sel(lat=slice(0, 90))).mean(dim='lat')
    sh_mean = da.sel(lat=slice(-90, 0)).weighted(weights.sel(lat=slice(-90, 0))).mean(dim='lat')
    
    # Create DataFrame for this species
    df = pd.DataFrame({
        'Year': decimal_years[year_mask],
        'Global': global_mean.values,
        'NH': nh_mean.values,
        'SH': sh_mean.values
    })
    
    # Store units for column headers
    units = da.attrs.get('units', 'Annual average (ppt)')
    all_data.append((species, df, units))

    ##########################################################################
    # Plotting
    ax = axs[i]
    ax.plot(global_mean['time'], global_mean, label='Global', color='black')
    ax.plot(nh_mean['time'], nh_mean, label='NH', color='red')
    ax.plot(sh_mean['time'], sh_mean, label='SH', color='blue')
    
    # Set title and labels
    ax.set_title(f"{species.upper()}", fontsize =9)
    ax.set_ylabel('Annual average (ppt)', fontsize=9)
    ax.grid(False)
    
    # Set x-axis to display integer years only
    ax.set_xticks(range(2000, 2023, 2))

# Add legend to first plot
axs[0].legend(loc='upper left', fontsize = 7)

# Set common xlabel
for ax in axs[4:]:
    ax.set_xlabel("Year",fontsize=9)

for ax in axs:
    ax.tick_params(axis='both', labelsize=7)  

plt.tight_layout()

# Save the figure
plt.savefig('FigureS2', dpi=300, bbox_inches='tight')
plt.show()
