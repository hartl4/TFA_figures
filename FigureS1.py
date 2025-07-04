#CODE TO PLOT FIGURE S1
##########################################################################################
#IMPORT MODULES
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.dates import date2num
import datetime
import subprocess
import numpy as np
##########################################################################################

# HCFC-123 DATA
agagedata = {
  "Years": [1996, 2000, 2003, 2004],
  "conc": [0.03, 0.05, 0.06, 0.064]
} 

hippodata_NH = {
  "Years": [2009, 2010, 2011],
  "conc": [0.182, 0.182, 0.182]
} #nh

hippodata_SH = {
  "Years": [2009, 2010, 2011],
  "conc": [0.1, 0.1, 0.1]
} #sh

NHincrease = { 
  'Years': [1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022], 
  'conc' : [0.0855, 0.091, 0.09646, 0.102248, 0.1092, 0.11648, 0.129584, 0.142688, 0.155792, 
    0.168896, 0.182, 0.182, 0.182,  0.182, 0.182, 0.182, 
   0.182, 0.182, 0.182, 0.182, 0.182, 0.182, 0.182, 0.182]} 

SHincrease = { 
  'Years': [1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022], 
  'conc' : [0.047, 0.05, 0.053, 0.05618, 0.06, 0.064, 0.0712, 0.0784, 0.0856, 
    0.0928, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1, 
   0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]} 

Global_av = { 
  'Years': [1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022],
  'conc' : [0.06627, 0.0705, 0.07473, 0.07921, 0.0846, 0.09024, 0.100392, 0.1105, 0.1207, 
    0.1308, 0.141, 0.141, 0.141,  0.141, 0.141, 0.141, 
   0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141, 0.141]} 

##############################################################################################################################################################
#HCFC-124 DATA
#LOAD AGAGE DATA (https://www-air.larc.nasa.gov/missions/agage/data/)

N_files = [f for f in os.listdir(directory_99) if f.startswith('N')] # Select Northern Hemisphere files
S_files = [f for f in os.listdir(directory_99) if f.startswith('S')] # Select Southern Hemisphere files

N_files_05 = [f for f in os.listdir(directory_05) if f.startswith('N')] # Select Northern Hemisphere files
S_files_05 = [f for f in os.listdir(directory_05) if f.startswith('S')] # Select Southern Hemisphere files

#CALCULATE HEMISPHERIC AND GLOBAL AVERAGES FOR  EACH YEAR
# Function to calculate annual averages for years 2000-2004
def calculate_annual_averages(files):
    total_per_year = {}
    count_per_year = {}
    
    # Iterate through each file
    for file_name in files:
        file_path = os.path.join(directory_99, file_name)
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.split()

                try:
                    # Extract the year, month, and monthly mean
                    year = int(columns[1])
                    mean = float(columns[13])              
                    
                    if not (mean != mean): # Skip NaN values
                        # Accumulate the total and count for the year
                        if year in total_per_year:
                            total_per_year[year] += mean
                            count_per_year[year] += 1
                        else:
                            total_per_year[year] = mean
                            count_per_year[year] = 1
                except (ValueError, IndexError):
                    continue

    # Calculate annual averages
    annual_averages = {}
    for year in sorted(total_per_year.keys()):
        if 1999< year <2005:
            total = total_per_year[year]
            count = count_per_year[year]
            average = total / count
            annual_averages[year] = average 
    return annual_averages

def calculate_annual_averages_05(files_05):
    total_per_year_05 = {}
    count_per_year_05 = {}
    
        # Iterate through each file
    for file_name_05 in files_05:
        file_path_05 = os.path.join(directory_05, file_name_05)
        with open(file_path_05, 'r') as file_05:
            # Iterate through each line of the file
            for line in file_05:
                # Split the line into columns based on whitespace
                columns_05 = line.split()

                try:
                    # Extract the year, month, and monthly mean
                    year_05 = int(columns_05[1])
                    mean_05 = float(columns_05[3])
                
                    # Skip NaN values
                    if not (mean_05 != mean_05):
                        # Accumulate the total and count for the year
                        if year_05 in total_per_year_05:
                            total_per_year_05[year_05] += mean_05
                            count_per_year_05[year_05] += 1
                        else:
                            total_per_year_05[year_05] = mean_05
                            count_per_year_05[year_05] = 1
                except (ValueError, IndexError):
                    # Skip lines that cannot be converted to float or do not have required columns
                    continue

    # Calculate annual averages for 2005-2022
    annual_averages_05 = {}
    #for year, total in total_per_year.items():
    for year_05 in sorted(total_per_year_05.keys()):
        if 2004< year_05:
            total_05 = total_per_year_05[year_05]
            count_05 = count_per_year_05[year_05]
            average_05 = total_05 / count_05
            annual_averages_05[year_05] = average_05
    
    return annual_averages_05

# Calculate Northern Hemisphere annual averages
N_annual_averages = calculate_annual_averages(N_files) 
N_annual_averages_05 = calculate_annual_averages_05(N_files_05)
N_sorted = sorted(N_annual_averages.items()) 
N05_sorted = sorted(N_annual_averages_05.items())
# Calculate Southern Hemisphere annual averages
S_annual_averages = calculate_annual_averages(S_files)
S_annual_averages_05 = calculate_annual_averages_05(S_files_05)
S_sorted = sorted(S_annual_averages.items())
S05_sorted = sorted(S_annual_averages_05.items())
#calculate global average
Global_annual_average = calculate_annual_averages(N_files + S_files)
Global_annual_average_05 = calculate_annual_averages_05(N_files_05 + S_files_05)

combined_NH = {**N_annual_averages, **N_annual_averages_05}
combined_SH = {**S_annual_averages, **S_annual_averages_05}
combined_G = {**Global_annual_average, **Global_annual_average_05}

#############################################################################################################
#HCFC-133A DATA
# Function to process hcfc-133a data from vollmer et al 
def process_species_data(vollmerdata):
    # Extract relevant columns
    years = vollmerdata.iloc[:, 0]
    NH_species = vollmerdata.iloc[:, 5]
    SH_species = vollmerdata.iloc[:, 6]

    # Create DataFrames with the year and concentration data
    df_NH = pd.DataFrame({'Year': years, 'Concentration': NH_species})
    df_SH = pd.DataFrame({'Year': years, 'Concentration': SH_species})
    df_NH['Year'] = df_NH['Year'].astype(str).str.extract(r'(\d{4})').astype(int)
    df_SH['Year'] = df_SH['Year'].astype(str).str.extract(r'(\d{4})').astype(int)

    # Group by Year and calculate annual averages (mean concentration)
    NH_annual_avg = df_NH.groupby('Year').mean().reset_index()
    SH_annual_avg = df_SH.groupby('Year').mean().reset_index()
    df_concat = pd.concat([NH_annual_avg, SH_annual_avg])

    # Filter for years greater than or equal to 2000
    df_concat = df_concat[df_concat['Year'] >= 2000]

    # Group by Year and calculate the mean concentration for each year
    global_averages = df_concat.groupby('Year').mean().reset_index()

    return years, NH_annual_avg, SH_annual_avg, global_averages

# Function to process hcfc-133a data from AGAGE data
def calculate_annual_averages_133(files, directory):
    total_per_year_133 = {}
    count_per_year_133 = {}

    for file_name in files:
        file_path = os.path.join(directory, file_name)
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.split()

                try:
                    # Extract the year and monthly mean
                    year_133 = int(columns[1])
                    mean_133 = float(columns[3])
                
                    if not (mean_133 != mean_133):  # Check if value is NaN
                        if year_133 in total_per_year_133:
                            total_per_year_133[year_133] += mean_133
                            count_per_year_133[year_133] += 1
                        else:
                            total_per_year_133[year_133] = mean_133
                            count_per_year_133[year_133] = 1
                except (ValueError, IndexError):
                    continue

    # Calculate annual averages
    annual_averages_133 = {}
    for year_133 in sorted(total_per_year_133.keys()):
        if year_133 > 2014:  # Use data from 2015 onwards
            total_133 = total_per_year_133[year_133]
            count_133 = count_per_year_133[year_133]
            annual_averages_133[year_133] = total_133 / count_133

    return annual_averages_133

#PROCESS VOLLMER DATA (2000-2014)
vollmer_file = '/pathto/HCFC133avollmer.xlsx'
vollmerdata = pd.read_excel(vollmer_file, sheet_name=0)
years, NH_annual_avg, SH_annual_avg, global_averages = process_species_data(vollmerdata)

#PROCESS AGAGE DATA
directory = '/pathto_agagedata'

N_files_133 = [f for f in os.listdir(directory) if f.endswith('a.txt')]  # Northern Hemisphere files
S_files_133 = [f for f in os.listdir(directory) if f.endswith('S.txt')]  # Southern Hemisphere files

N_annual_averages_133 = calculate_annual_averages_133(N_files_133, directory)  # NH averages from 2015
S_annual_averages_133 = calculate_annual_averages_133(S_files_133, directory)  # SH averages from 2015
Global_annual_average_133 = calculate_annual_averages_133(N_files_133 + S_files_133, directory)

# Combine data up to 2014 from Vollmer and from 2015 onwards from AGAGE
combined_NH_133 = pd.concat([
    pd.DataFrame({'Year': NH_annual_avg['Year'], 'Concentration': NH_annual_avg['Concentration']}),
    pd.DataFrame(list(N_annual_averages_133.items()), columns=['Year', 'Concentration']) 
])

combined_SH_133 = pd.concat([pd.DataFrame({'Year': SH_annual_avg['Year'], 'Concentration': SH_annual_avg['Concentration']}),
                    pd.DataFrame(list(S_annual_averages_133.items()), columns=['Year', 'Concentration'])])

combined_global_133 = pd.concat([pd.DataFrame({'Year': global_averages['Year'], 'Concentration': global_averages['Concentration']}),pd.DataFrame(list(Global_annual_average_133.items()), columns=['Year', 'Concentration'])])

################################################################################################################
#PLOTTING
################################################################################################################
# Create a figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(9, 3))  # 3 rows, 1 column of subplots

# Plot HCFC-123 
axes[0].plot(SHincrease["Years"], SHincrease["conc"], linestyle='-', color='blue', label='SH')
axes[0].plot(NHincrease["Years"], NHincrease["conc"], linestyle='-', color='red', label='NH')
axes[0].plot(Global_av["Years"], Global_av["conc"], linestyle='-', color='black', label='Global')
axes[0].set_title('HCFC-123',fontsize=9)
axes[0].set_xlabel('Year',fontsize=9)
axes[0].set_ylabel('Annual Average (ppt)',fontsize=9)
axes[0].tick_params(axis='both', labelsize=7)  

# Plot HCFC-124 
axes[1].plot(combined_NH.keys(), combined_NH.values(), linestyle='-', label='HCFC-124 NH', color = 'red')
axes[1].plot(combined_SH.keys(), combined_SH.values(), linestyle='-', label='HCFC-124 SH', color = 'blue')
axes[1].plot(combined_G.keys(), combined_G.values(), linestyle='-', label=' HCFC-124Global', color = 'black')
axes[1].set_title('HCFC-124',fontsize=9)
axes[1].set_xlabel('Year',fontsize=9)
axes[1].tick_params(axis='both', labelsize=7) 
axes[1].set_ylabel('Annual Average (ppt)',fontsize=9)

# Plot HCFC-133a 
axes[2].plot(combined_NH_133['Year'], combined_NH_133['Concentration'], linestyle='-', label='HCFC-133a NH', color='red')
axes[2].plot(combined_SH_133['Year'], combined_SH_133['Concentration'], linestyle='-', label='HCFC-133a SH', color='blue')
axes[2].plot(combined_global_133['Year'], combined_global_133['Concentration'], linestyle='-', label='HCFC-133a Global', color='black')
axes[2].set_title('HCFC-133a',fontsize=9)
axes[2].set_xlabel('Year',fontsize=9)
axes[2].set_ylabel('Annual Average (ppt)',fontsize=9)
axes[2].tick_params(axis='both', labelsize=7)  

# Create horizontal legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=7, frameon=False,bbox_to_anchor=(0.54, -0.03))

plt.tight_layout()
plt.savefig('FigureS1')
plt.show()