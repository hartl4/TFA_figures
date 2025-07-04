# Description: This code uses measured data from Pickard et al, 2020, Hartz et al. 2023, and Freeling et al. 2020 
                # as well as TFA fluxes (kg/gricbox/month) at the locations of the measurements.
# Author: Lucy Hart
# Date: 09/06/25

#IMPORT MODULES
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from mpl_toolkits.basemap import shiftgrid
import cartopy.crs as ccrs
import csv
from matplotlib.ticker import MaxNLocator
import cmcrameri.cm as cmc

#######################################################################################
#PREP MEASURED ICE CORE DATA FOR COMPARISON TO MODEL
#######################################################################################

#CALCULATE 5 YEAR ROLLING AVERAGE FROM PICKARD ET AL DATA
#file_path = '/filepathto/Pickarddata2.csv'
file_path = '/home/hartl4/luna/Cryo-Chem/hartl4/graphs/hcfcs/Pickarddata2.csv'

#define indexes of ice core locations from files at T63 resolution
locations = [
    {'name': 'oxford', 'lat_index': 91, 'lon_index': 57},
    {'name': 'devon', 'lat_index': 88, 'lon_index': 52},
    {'name': 'svalbard', 'lat_index': 89, 'lon_index': 105}] 
    
# Define window size for moving average
window_size = 5
moving_averages_location = {loc['name']: [] for loc in locations}
years_moving = []
df = pd.read_csv(file_path)
csv_years = df['year'].tolist()
devon = df['Devon'].tolist()
oxford = df['Oxford'].tolist()
svalbard = df['Svalbard'].tolist()
location_data = {
    'svalbard': svalbard,
    'devon': devon,
    'oxford': oxford
}

# Loop through the array and calculate moving average
for i in range(len(csv_years) - window_size + 1):
    for location in locations:
        name = location['name']
        window = location_data[name][i: i + window_size]
        window_average = np.mean(window)
        moving_averages_location[name].append(window_average)   
    years_moving.append(csv_years[i + window_size // 2])
    
# Function to calculate flux mean and standard deviation for a central lat/lon and its surrounding 9-box region
def calculate_flux(lat_index, lon_index, TFA_wetdep, TFA_drydep, surfarea):
    # Sum the arrays over the time and lev dimensions to calculate annual flux
    tfa_wet = TFA_wetdep.sum(dim=('time', 'lev')) #calc wet depoosition flux in kg/gridbox 
    tfa_dry = TFA_drydep.sum(dim=('time', 'lev')) #calc dry depoosition flux in kg/gridbox 
    surf_time = surfarea.mean(dim='time').squeeze()

    # Define the indices for the 3x3 grid around the central lat/lon
    lat_indices = range(lat_index - 1, lat_index + 2)  # Get lat indices of gridbox and surrounding boxes for S.D. calc
    lon_indices = range(lon_index - 1, lon_index + 2)  # Get lON indices of gridbox and surrounding boxes for S.D. calc
    tfa_wet_flux = tfa_wet.isel(lat=lat_indices, lon=lon_indices).values
    tfa_dry_flux = tfa_dry.isel(lat=lat_indices, lon=lon_indices).values

    # Calculate gridbox area
    area = surf_time.isel(lat=lat_indices, lon=lon_indices).values 

    # Calculate the total TFA deposition flux
    tfa_total_dep = abs(tfa_dry_flux) + abs(tfa_wet_flux)
    tfa_dep = (tfa_total_dep / area) * 1e12  # Convert to ng/m²/yr

    # Calculate mean and standard deviation across the 3x3 grid
    tfa_mean = np.mean(tfa_dep)  
    tfa_std = np.std(tfa_dep)    

    return {
        'tfa_mean': tfa_mean,
        'tfa_std': tfa_std
    }

# Dictionary to store flux data by year for each location
flux_by_location = {loc['name']: {'mean': [], 'std': []} for loc in locations}
years = []

####################################################################################
#CALCULATE MODEL DEPOSITION FLUXES AT LOCATIONS OF ICE CORES (MEAN +- S.D.)
####################################################################################

#files =  [f"/path to/tend4d_{year}_latlon.nc" for year in range(2000, 2020)]
#files_surf =  [f"/pathto/surf_{year}_latlon.nc" for year in range(2000, 2020)]
#file_path = '/path to/Pickarddata2.csv'
files =  [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/cmiplus/base/tend4d_{year}_latlon.nc" for year in range(2000, 2020)]
files_surf =  [f"/home/hartl4/luna/Cryo-Chem/hartl4/outputs/chpt1/simple/surf_{year}_latlon.nc" for year in range(2000, 2020)]
file_path = '/home/hartl4/luna/Cryo-Chem/hartl4/graphs/hcfcs/Pickarddata2.csv'
# Loop through the files for each year
for file, file_surf in zip(files, files_surf):
    year = file.split('4d_')[1].split('_lat')[0]
    years.append(int(year))

    # Open the dataset
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    sa = xr.open_dataset(file_surf, engine='netcdf4', decode_times=False)

    # Extract the relevant variables from the dataset
    surfarea = sa['area']  
    TFA_wetdep = ds['tfa_wetcnv'] #kg/gridbox/month
    TFA_drydep = ds['tfa_drydep'] #kg/gridbox/month

    # Loop through each location and calculate flux
    for location in locations:
        lat_index = location['lat_index']
        lon_index = location['lon_index']
        flux_data = calculate_flux(lat_index, lon_index, TFA_wetdep, TFA_drydep, surfarea)

        # Append the mean and std flux data for the location
        flux_by_location[location['name']]['mean'].append(flux_data['tfa_mean'])
        flux_by_location[location['name']]['std'].append(flux_data['tfa_std'])

#read pickard et al. data from .csv
csv_years = df['year'].tolist()
devon_pickard = df['Devon'].tolist()
oxford_pickard = df['Oxford'].tolist()
svalbard_ice = df['Svalbard'].tolist()

# Filter the CSV data to include only the years from 2000 onwards
csv_filtered_years = [year for year in csv_years if 2000 <= year <= 2023]
devon_pickard_filtered = [devon_pickard[i] for i, year in enumerate(csv_years) if 2000<= year <= 2023]
oxford_pickard_filtered = [oxford_pickard[i] for i, year in enumerate(csv_years) if 2000<= year <= 2023]
svalbard_ice_filtered = [svalbard_ice[i] for i, year in enumerate(csv_years) if 2000<= year <= 2023]

########################################################################
#prep the FREELING et al rainwater concentration data
########################################################################

# Load the data
#file = 'pathto/6_Einzel_Misch_worst_case_mass.xlsx'
#file2 = '/pathto/5_Einzel_Misch_best_case_mass.xlsx'
file = '/home/hartl4/luna/Cryo-Chem/hartl4/graphs/freeling/6_Einzel_Misch_worst_case_mass.xlsx'
file2 = '/home/hartl4/luna/Cryo-Chem/hartl4/graphs/freeling/5_Einzel_Misch_best_case_mass.xlsx'
df = pd.read_excel(file)
df2 = pd.read_excel(file2)
frames = [df,df2]
dfs = pd.concat(frames)

# Reshape the DataFrame
df_melted = dfs.melt(id_vars=["Station"], var_name="year_month", value_name="value").dropna()
df_melted["year_month"] = pd.to_datetime(df_melted["year_month"], format="%Y_%m")

# Filter the data for Feb 2018 to Jan 2019
start_date = "2018-02"
end_date = "2019-01"
df_filtered = df_melted[(df_melted["year_month"] >= start_date) & (df_melted["year_month"] <= end_date)]

#####################################################################################
#CALCULATE MODEL RAINWATER CONCENTRATIONS FOR COMPARISON TO MEASUREMENT
#####################################################################################
'''
# Sort by date to ensure correct order
df_filtered = df_filtered.sort_values("year_month")
fn ='/path_to/tend4d_2018_latlon.nc'
ds = xr.open_dataset(fn,engine='netcdf4',decode_times=False)
fn2 ='/path_to/tend4d_2018_latlon.nc'
ds2 = xr.open_dataset(fn2,engine='netcdf4',decode_times=False)
path ='/path_to/surf_2018_latlon.nc'
sa = xr.open_dataset(path,engine='netcdf4',decode_times=False)
path2 ='/path_to/surf_2018_latlon.nc'
sa2 = xr.open_dataset(path2,engine='netcdf4',decode_times=False)
file = '/path_to/7_Misch_Misch_best_case_mass.xlsx'
df = pd.read_excel(file)
'''
# Sort by date to ensure correct order
df_filtered = df_filtered.sort_values("year_month")
#fn ='/home/hartl4/luna/Cryo-Chem/hartl4/outputs/hfo/tend4d_2018_latlon.nc'
fn ='/home/hartl4/luna/Cryo-Chem/hartl4/outputs/expts/hfo_high/tend4d_2018_latlon.nc'
ds = xr.open_dataset(fn,engine='netcdf4',decode_times=False)
fn2 ='/home/hartl4/luna/Cryo-Chem/hartl4/outputs/expts/base/tend4d_2018_latlon.nc'
ds2 = xr.open_dataset(fn2,engine='netcdf4',decode_times=False)
path ='/home/hartl4/luna/Cryo-Chem/hartl4/outputs/chpt1/simple/surf_2018_latlon.nc'
sa = xr.open_dataset(path,engine='netcdf4',decode_times=False)
path2 ='/home/hartl4/luna/Cryo-Chem/hartl4/outputs/chpt1/simple/surf_2018_latlon.nc'
sa2 = xr.open_dataset(path2,engine='netcdf4',decode_times=False)
file = '/home/hartl4/luna/Cryo-Chem/hartl4/graphs/freeling/7_Misch_Misch_best_case_mass.xlsx'
df = pd.read_excel(file)

#dictionary to store lat and lon of measurement sites
sites = [
    {'name': 'SW', 'lat': 54.5275, 'lon': 9.5487},
    {'name': 'GW', 'lat': 54.0967, 'lon': 13.4056},
    {'name': 'PD', 'lat': 52.3813, 'lon': 13.0622},
    {'name': 'BR', 'lat': 51.7986, 'lon': 10.6183},
    {'name': 'ES', 'lat': 51.4041, 'lon': 6.9677},
    {'name': 'WK', 'lat': 50.4973	, 'lon': 9.9428},
    {'name': 'SU', 'lat': 48.8281, 'lon': 9.2},
    {'name': 'MO', 'lat': 48.244, 'lon': 11.553	}]

df_melted = df.melt(id_vars=["Station"], var_name="year_month", value_name="value")
df_melted["year_month"] = pd.to_datetime(df_melted["year_month"], format="%Y_%m") # convert year_monthg to datetime

# Filter the data for Feb 2018 to Jan 2019
start_date = "2018-02"
end_date = "2019-01"
df_filtered = df_melted[(df_melted["year_month"] >= start_date) & (df_melted["year_month"] <= end_date)]
df_filtered = df_filtered.sort_values("year_month")

# CALCULATE TFA CONCENTRATION AT EACH SITE PER MONTH AND ADD TO DICTIONARY
# extract variables
tfa_wet = ds2['tfa_wetcnv'].sum(dim='lev')
tfa_dry = ds2['tfa_drydep'].sum(dim='lev')
area = sa2['area'].mean(dim='lev')
monthly_values = {}

tfa_wet_hfo = ds['tfa_wetcnv'].sum(dim='lev')
tfa_dry_hfo = ds['tfa_drydep'].sum(dim='lev')
area_hfo = sa['area'].mean(dim='lev')
monthly_values_hfo = {}

# Convert hours to dates and then to months
reference_date = np.datetime64('2018-02-01T01:00:00')
dates = [reference_date + np.timedelta64(int(hours), 'h') for hours in ds.time.values]
months = [date.astype('datetime64[M]').astype(int) % 12 + 1 for date in dates]

# Create a mapping from month number to time indices
month_to_idx = {}
for idx, month in enumerate(months):
    if month not in month_to_idx:
        month_to_idx[month] = []
    month_to_idx[month].append(idx)

# Extract values for each month and site
for month in range(1, 13):
    if month not in month_to_idx:
        print(f"Warning: Month {month} not found in dataset")
        monthly_values[month] = [float('nan')] * len(sites)
        continue

    # Use the first occurrence of the month if there are multiple
    time_idx = month_to_idx[month][0]
    site_values = []
    site_values_hfo = []

#extract site values for base case
    for site in sites:
        try:
            # Find nearest grid points for site coordinates
            lat = site['lat']
            lon = site['lon']

            if site == 'BR':
                tfa_data_monthly = tfa_wet.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest') + tfa_dry.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
                area_data_monthly = area.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
                print(site)
            else:
            # Select data for current month and site location separately for tfa and prec before doing division
                tfa_data_monthly = tfa_wet.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
                area_data_monthly = area.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
            
            # Compute concentration and handle potential division issues
            with np.errstate(divide='ignore', invalid='ignore'):
                month_data = np.divide(tfa_data_monthly*1E9, area_data_monthly)
            
            # Flatten and take absolute value
            value = float(abs(month_data.values.flatten()[0]))
            site_values.append(value)

        except (KeyError, IndexError) as e:
            print(f"Warning: Error processing site {site} for month {month}: {str(e)}")
            site_values.append(float('nan'))

    monthly_values[month] = site_values

#extract site values for hfo run
    for site in sites:
        try:
            # Find nearest grid points for site coordinates
            lat = site['lat']
            lon = site['lon']

            if site == 'BR':
                tfa_data_monthly_hfo = tfa_wet_hfo.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest') + tfa_dry_hfo.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
                area_data_monthly_hfo = area_hfo.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
            else:
            # Select data for current month and site location separately for tfa and prec before doing division
                tfa_data_monthly_hfo = tfa_wet_hfo.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
                area_data_monthly_hfo = area_hfo.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
            
            # Compute concentration and handle potential division issues
            with np.errstate(divide='ignore', invalid='ignore'):
                month_data_hfo = np.divide(tfa_data_monthly_hfo*1E9, area_data_monthly_hfo)
            
            # Flatten and take absolute value
            value = float(abs(month_data_hfo.values.flatten()[0]))
            site_values_hfo.append(value)

        except (KeyError, IndexError) as e:
            print(f"Warning: Error processing site {site} for month {month}: {str(e)}")
            site_values_hfo.append(float('nan'))

    monthly_values_hfo[month] = site_values_hfo

# Define the ordered months we want to plot (Feb 2018 - Jan 2019)
ordered_months = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]

# Prepare Freeling data for plotting
# Group the filtered data by month
freeling_groups = df_filtered.groupby(df_filtered['year_month'].dt.month)

# Prepare lists to store Freeling data for each month
freeling_monthly_values = {}
for month, group in freeling_groups:
    freeling_monthly_values[month] = group['value'].tolist()

# Filter Nans from feb data
if 2 in freeling_monthly_values:
    feb_data = freeling_monthly_values[2]
    freeling_monthly_values[2] = [val for val in feb_data if not (isinstance(val, (int, float)) and np.isnan(val))]

######################################################################################
#PLOTTING
######################################################################################

# Create figure and axis
fig, ax = plt.subplots(figsize=(15, 7))

plot_data = []
positions = []
labels = []
colors = []

# Process each month in our desired order
for i, month in enumerate(ordered_months):
    # Add model data
    if month in monthly_values:
        # Filter out NaN values from model data
        model_data = [val for val in monthly_values[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if model_data:
            plot_data.append(model_data)
            positions.append(i * 3)
            labels.append(f'{month} (Model)')
            colors.append('lightblue')
            print(f"Added model data for month {month}: {len(model_data)} values")
        else:
            print(f"WARNING: Only NaN values in model data for month {month}")
    else:
        print(f"WARNING: No model data for month {month}")
    
    # Add Freeling measured data
    if month in freeling_monthly_values:
        # Filter out NaN values from Freeling data
        freeling_data = [val for val in freeling_monthly_values[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if freeling_data:  # Only add if there's non-NaN data
            plot_data.append(freeling_data)
            positions.append(i * 3 + 1)
            labels.append(f'{month} (Freeling)')
            colors.append('lightgreen')
            print(f"Added Freeling data for month {month}: {len(freeling_data)} values")
        else:
            print(f"WARNING: Only NaN values in Freeling data for month {month}")
    else:
        print(f"WARNING: No Freeling data for month {month}")

# Create the subplots with 2x2 grid
fig, axs = plt.subplots(2, 2, figsize=(16, 8), constrained_layout=True)  # Wider figure for rectangular plots

# Define labels for each location
location_labels = {
    'devon': 'Devon Ice Cap',
    'oxford': 'Mt. Oxford Icefield',
    'svalbard': 'Svalbard'
}

# Left-hand panels (Devon and Oxford data)
# Devon Ice Cap (top left)
min_year = int(min(min(csv_filtered_years), min(years_moving), min(years)))
max_year = int(max(max(csv_filtered_years), max(years_moving), max(years))) + 1
whole_years = np.arange(min_year, max_year, 2)  # Step of 5 years 
axs[0, 0].set_xticks(np.arange(min_year, max_year + 1, 2))  # every 2 years
axs[0, 0].set_xticklabels([str(year) for year in np.arange(min_year, max_year + 1, 2)])
axs[0, 0].plot(csv_filtered_years, devon_pickard_filtered, label='Measured', color='black', linestyle='--')
axs[0, 0].plot(years_moving, moving_averages_location['devon'], label='Measured (5 Year Moving Average)', color='black', linestyle='-')
axs[0, 0].errorbar(years, flux_by_location['devon']['mean'], yerr=flux_by_location['devon']['std'], 
                   label='Model: BASE (Mean ± S.D.)', marker='o', color='orange', capsize=3)
axs[0, 0].set_title('(a) Devon Ice Cap', fontsize=14)
axs[0, 0].set_ylabel('TFA Deposition (ng/m²/yr)', fontsize=14)
axs[0, 0].legend(loc='upper left')

# Mt. Oxford Icefield (bottom left)
axs[1, 0].set_xticks(np.arange(min_year, max_year + 1, 2))  # every 2 years
axs[1, 0].set_xticklabels([str(year) for year in np.arange(min_year, max_year + 1, 2)])
axs[1, 0].plot(csv_filtered_years, oxford_pickard_filtered, label='Measured', color='black', linestyle='--')
axs[1, 0].plot(years_moving, moving_averages_location['oxford'], label='Measured (5 Year Moving Average)', color='black', linestyle='-')
axs[1, 0].errorbar(years, flux_by_location['oxford']['mean'], yerr=flux_by_location['oxford']['std'], 
                   label='Model: BASE (Mean ± S.D.)', marker='o', color='green', capsize=3)
axs[1, 0].set_title('(b) Mt. Oxford Icefield', fontsize=14)
axs[1, 0].set_ylabel('TFA Deposition (ng/m²/yr)', fontsize=14)
axs[1, 0].set_xlabel('Year', fontsize=12)
axs[1, 0].legend(loc='upper left')
axs[0, 1].set_xticks(np.arange(min_year, max_year + 1, 2))  # every 2 years
axs[0, 1].set_xticklabels([str(year) for year in np.arange(min_year, max_year + 1, 2)])

# Bottom-right panel: Svalbard data
axs[0, 1].plot(csv_filtered_years, svalbard_ice_filtered, label='Measured', color='black', linestyle='--')
axs[0, 1].plot(years_moving, moving_averages_location['svalbard'], label='Measured (5 Year Moving Average)', color='black', linestyle='-')
axs[0, 1].errorbar(years, flux_by_location['svalbard']['mean'], yerr=flux_by_location['svalbard']['std'], 
                   label='Model: BASE (Mean ± S.D.)', marker='o', color='blue', capsize=3)
axs[0, 1].set_title('(c) Lomonosovfonna', fontsize=14)
axs[0, 1].set_ylabel('TFA Deposition (ng/m²/yr)', fontsize=14)
axs[0, 1].set_xlabel('Year', fontsize=14)
axs[0, 1].legend(loc='upper left')

#TOP RIGHT PANEL (FREELING FLUX)
# Define the ordered months to plot (Feb 2018 - Jan 2019)
ordered_months = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]
# Group the filtered data by month
freeling_groups = df_filtered.groupby(df_filtered['year_month'].dt.month)
# Define the ordered months to plot (Feb 2018 - Jan 2019)
ordered_months = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]

# Prepare lists to store Freeling data for each month
freeling_monthly_values = {}
for month, group in freeling_groups:
    freeling_monthly_values[month] = group['value'].tolist()

plot_data = []
positions = []
labels = []
colors = []
box_colors = {'Model': 'lightblue', 'Model HFO': 'darkblue', 'Freeling': 'grey'}

# Process each month in our desired order
for i, month in enumerate(ordered_months):
    # Add original model data
    if month in monthly_values:
        # Filter out NaN values from model data
        model_data = [val for val in monthly_values[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if model_data: 
            plot_data.append(model_data)
            positions.append(i * 4)
            labels.append(f'{month} (Model)')
            colors.append(box_colors['Model'])
            print(f"Added model data for month {month}: {len(model_data)} values")
        else:
            print(f"WARNING: Only NaN values in model data for month {month}")
    else:
        print(f"WARNING: No model data for month {month}")
    
    # Add HFO model data
    if month in monthly_values_hfo:
        # Filter out NaN values from HFO model data
        model_data_hfo = [val for val in monthly_values_hfo[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if model_data_hfo:  # Only add if there's non-NaN data
            plot_data.append(model_data_hfo)
            positions.append(i * 4 + 1)
            labels.append(f'{month} (Model HFO)')
            colors.append(box_colors['Model HFO'])
            print(f"Added HFO model data for month {month}: {len(model_data_hfo)} values")
        else:
            print(f"WARNING: Only NaN values in HFO model data for month {month}")
    else:
        print(f"WARNING: No HFO model data for month {month}")
    
    # Add Freeling data
    if month in freeling_monthly_values:
        # Filter out NaN values from Freeling data
        freeling_data = [val for val in freeling_monthly_values[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if freeling_data:  # Only add if there's non-NaN data
            plot_data.append(freeling_data)
            positions.append(i * 4 + 2)
            labels.append(f'{month} (Freeling)')
            colors.append(box_colors['Freeling'])
            print(f"Added Freeling data for month {month}: {len(freeling_data)} values")
        else:
            print(f"WARNING: Only NaN values in Freeling data for month {month}")
    else:
        print(f"WARNING: No Freeling data for month {month}")
        
# Add a border around the axis
from matplotlib.patches import Rectangle
rect = Rectangle((0, 0), 1, 1, transform=axs[1,1].transAxes, edgecolor='black', facecolor='none', lw=1)
axs[1,1].add_patch(rect)

# Create boxplot
bp = axs[1,1].boxplot(plot_data,
                positions=positions,
                patch_artist=True,
                medianprops=dict(color="black", linewidth=1.5),
                flierprops=dict(marker='o', markerfacecolor='gray', markersize=4))

# Set colors for boxes
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# Customize plot
axs[1,1].set_title('(d) Monthly TFA Flux: Model vs Measurements (Feb 2018 - Jan 2019)',fontsize=14)
axs[1,1].set_ylabel('TFA Flux ((\u03bcg/m\u00b2)',fontsize=14)
axs[1,1].set_yscale('log')

# Set x-ticks at the middle of each month's group of boxes
month_positions = [i * 4 + 1 for i in range(len(ordered_months))]
axs[1,1].set_xticks(month_positions)

# Use proper month names in the correct order
month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
month_labels = [month_names[m-1] for m in ordered_months]  
axs[1,1].set_xticklabels(month_labels)

# Add vertical dotted lines between months
for i in range(1, len(ordered_months)):
    x_pos = i * 4 - 0.5  # Position between months
    axs[1,1].axvline(x=x_pos, color='gray', linestyle=':', alpha=0.7)

# Add grid
axs[1,1].grid(True, axis='y', linestyle='--', alpha=0.7)

# Create a custom legend with colored boxes
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=box_colors['Model'], edgecolor='black', label='Model (run BASE)'),
    Patch(facecolor=box_colors['Model HFO'], edgecolor='black', label='Model (run HFO high)'),
    Patch(facecolor=box_colors['Freeling'], edgecolor='black', label='Measured (Freeling et al. (2020))')
]
    
# Add the custom legend
axs[1,1].legend(handles=legend_elements, loc='upper right', 
          bbox_to_anchor=(1, 1), framealpha=0.5)

# Save and show the plot
plt.savefig('Figure3.png', dpi=300, bbox_inches = 'tight')
plt.show()
