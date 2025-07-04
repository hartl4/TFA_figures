##########################################################################################
#CODE TO PLOT FIGURE S7
##########################################################################################
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
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

########################################################################
#prep the FREELING et al rainwater concentration data
########################################################################

# Load the Freeling et al data
#file = 'pathto/6_Einzel_Misch_worst_case_mass.xlsx'
#file2 = '/pathto/5_Einzel_Misch_best_case_mass.xlsx'
file = 'freelingdata.xlsx'
file2 = 'freelingdata.xlsx'
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
#CALCULATE MODEL RAINWATER CONCENTRATIONS FROM MODEL FOR COMPARISON TO MEASUREMENT
#####################################################################################

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
df_melted["year_month"] = pd.to_datetime(df_melted["year_month"], format="%Y_%m") # convert year_mont to datetime

# Filter the data for Feb 2018 to Jan 2019
start_date = "2018-02"
end_date = "2019-01"
df_filtered = df_melted[(df_melted["year_month"] >= start_date) & (df_melted["year_month"] <= end_date)]
df_filtered = df_filtered.sort_values("year_month")

# CALCULATE TFA CONCENTRATION AT EACH SITE PER MONTH AND ADD TO DICTIONARY
tfa_wet = ds2['tfa_wetcnv'].sum(dim='lev') 
tfa_dry = ds2['tfa_drydep'].sum(dim='lev')
area = sa2['area'].mean(dim='lev') #m2
monthly_values = {}

tfa_wet_hfo = ds['tfa_wetcnv'].sum(dim='lev') # wet deposition
tfa_dry_hfo = ds['tfa_drydep'].sum(dim='lev') # dry deposition
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
            # Select data for current month and site location 
                tfa_data_monthly = tfa_wet.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
                area_data_monthly = area.isel(time=time_idx).sel(lat=lat, lon=lon, method='nearest')
            
            # Calculate concentration as TFA mass/area
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
freeling_groups = df_filtered.groupby(df_filtered['year_month'].dt.month)

# Prepare lists to store Freeling data for each month
freeling_monthly_values = {}
for month, group in freeling_groups:
    freeling_monthly_values[month] = group['value'].tolist()

# Filter Nans from feb data
if 2 in freeling_monthly_values:
    feb_data = freeling_monthly_values[2]
    freeling_monthly_values[2] = [val for val in feb_data if not (isinstance(val, (int, float)) and np.isnan(val))]

#####################################################################################
#PLOTTING
#####################################################################################
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

# PLOTTING - SINGLE PLOT

fig, ax = plt.subplots(figsize=(12, 6))  # Just one plot now

# Group Freeling data by month
freeling_groups = df_filtered.groupby(df_filtered['year_month'].dt.month)

# Define the ordered months (Feb 2018 - Jan 2019)
ordered_months = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1]

# Prepare lists
freeling_monthly_values = {}
for month, group in freeling_groups:
    freeling_monthly_values[month] = group['value'].tolist()

plot_data = []
positions = []
labels = []
colors = []
box_colors = {'Model': 'lightblue', 'Model HFO': 'darkblue', 'Freeling': 'grey'}

# Populate boxplot data
for i, month in enumerate(ordered_months):
    if month in monthly_values:
        model_data = [val for val in monthly_values[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if model_data:
            plot_data.append(model_data)
            positions.append(i * 4)
            labels.append(f'{month} (Model)')
            colors.append(box_colors['Model'])

    if month in monthly_values_hfo:
        model_data_hfo = [val for val in monthly_values_hfo[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if model_data_hfo:
            plot_data.append(model_data_hfo)
            positions.append(i * 4 + 1)
            labels.append(f'{month} (Model HFO)')
            colors.append(box_colors['Model HFO'])

    if month in freeling_monthly_values:
        freeling_data = [val for val in freeling_monthly_values[month] if not (isinstance(val, (int, float)) and np.isnan(val))]
        if freeling_data:
            plot_data.append(freeling_data)
            positions.append(i * 4 + 2)
            labels.append(f'{month} (Freeling)')
            colors.append(box_colors['Freeling'])

# Add border
rect = Rectangle((0, 0), 1, 1, transform=ax.transAxes, edgecolor='black', facecolor='none', lw=1)
ax.add_patch(rect)

# Boxplot
bp = ax.boxplot(plot_data,
                positions=positions,
                patch_artist=True,
                medianprops=dict(color="black", linewidth=1.5),
                flierprops=dict(marker='o', markerfacecolor='gray', markersize=4))

# Set colors
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# Customize
ax.set_title('Monthly TFA Flux: Model vs Measurements (Feb 2018 - Jan 2019)', pad=20, fontsize = 10)
ax.set_ylabel('TFA Flux (μg/m²)', fontsize = 10)
ax.set_yscale('log')

# X-ticks and labels
month_positions = [i * 4 + 1 for i in range(len(ordered_months))]
month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
month_labels = [month_names[m-1] for m in ordered_months]
ax.set_xticks(month_positions)
ax.set_xlabel('Month', fontsize = 10)
ax.set_xticklabels(month_labels)

# Add vertical lines
for i in range(1, len(ordered_months)):
    x_pos = i * 4 - 0.5
    ax.axvline(x=x_pos, color='gray', linestyle=':', alpha=0.7)

# Grid
ax.grid(True, axis='y', linestyle='--', alpha=0.7)

# Legend
legend_elements = [
    Patch(facecolor=box_colors['Model'], edgecolor='black', label='Model (run BASE)'),
    Patch(facecolor=box_colors['Model HFO'], edgecolor='black', label='Model (run HFO high)'),
    Patch(facecolor=box_colors['Freeling'], edgecolor='black', label='Measured (Freeling et al. (2020))')
]
ax.legend(handles=legend_elements, loc='upper right', framealpha=0.5)

# Save and show
plt.savefig('Figures7.png', dpi=300, bbox_inches='tight')
plt.show()
