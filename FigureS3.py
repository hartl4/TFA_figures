################################################################################################
#CODE TO PLOT FIGURE S3
################################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

################################################################################################
#READ ANAEST DATA FROM VOLLMER ET ALFOR PLOTTING 
# Function to process the data for anesthetic species
def process_species_data(xl_data, species_list, species_names):
    processed_data = {}
    years = xl_data['t'].astype(int)

    for species, name in zip(species_list, species_names):
        # Create DataFrame with the year and concentration data
        df_mean = pd.DataFrame({'Year': years, 'Concentration': species})

        # Filter for years >= 1999 and < 2015
        df_mean = df_mean[(df_mean['Year'] >= 1999) & (df_mean['Year'] < 2015)]

        # Clean concentration values
        df_mean['Concentration'] = df_mean['Concentration'].replace('	  NaN', 0)
        df_mean['Concentration'] = df_mean['Concentration'].replace(r'\s*NaN\s*', 0, regex=True)
        df_mean['Concentration'] = pd.to_numeric(df_mean['Concentration'], errors='coerce').fillna(0)

        # Hold 2014 value constant for 2015–2022
        value_2014 = df_mean.loc[df_mean['Year'] == 2014, 'Concentration'].mean()

        extension_years = pd.DataFrame({
            'Year': range(2015, 2023),
            'Concentration': value_2014
        })
        df_extended = pd.concat([df_mean, extension_years], ignore_index=True)

        # Group by Year and calculate mean
        global_averages = df_extended.groupby('Year').mean().reset_index()

        # Store data
        processed_data[name] = {
            'years': df_extended['Year'],
            'concentrations': df_extended['Concentration'],
            'global_averages': global_averages,
        }

    return processed_data

################################################################################################
# Load the data
file = 'pathto/anaesth_vollmer.xlsx'
xl_data = pd.read_excel(file, skiprows=9)

# Extract species data
halothane = xl_data['HAmfBc']
Halothane = halothane*1e-3 #convert units to ppt
Isoflurane = xl_data['IFmfBc']
Desflurane = xl_data['DFmfBc']
Sevoflurane = xl_data['SFmfBc']

# List of species and their names
species_list = [Halothane, Isoflurane, Desflurane, Sevoflurane]
species_names = ['Halothane', 'Isoflurane', 'Desflurane', 'Sevoflurane']

# Process species data
species_data = process_species_data(xl_data, species_list, species_names)

################################################################################################
# Create a figure with 2x2 subplots (one for each species)
fig, axes = plt.subplots(1, 4, figsize=(12, 3)) 

# Plot Halothane
axes[0].plot(species_data['Halothane']['years'], species_data['Halothane']['concentrations'],color='black', linestyle='-')
axes[0].set_title('Halothane',fontsize=9)
axes[0].set_xlabel('Year',fontsize=9)
axes[0].set_ylabel('Annual Average (ppt)',fontsize=9)
axes[0].tick_params(axis='both', labelsize=7) 

# Plot Isoflurane
axes[ 1].plot(species_data['Isoflurane']['years'], species_data['Isoflurane']['concentrations'], color='black', linestyle='-')
axes[ 1].set_title('Isoflurane',fontsize=9)
axes[ 1].set_xlabel('Year',fontsize=9)
axes[ 1].set_ylabel('Annual Average (ppt)',fontsize=9)
axes[1].tick_params(axis='both', labelsize=7) 

# Plot Desflurane
axes[ 2].plot(species_data['Desflurane']['years'], species_data['Desflurane']['concentrations'], color='black',linestyle='-')
axes[ 2].set_title('Desflurane',fontsize=9)
axes[ 2].set_xlabel('Year',fontsize=9)
axes[ 2].set_ylabel('Annual Average (ppt)',fontsize=9)
axes[2].tick_params(axis='both', labelsize=7) 

# Plot Sevoflurane
axes[3].plot(species_data['Sevoflurane']['years'], species_data['Sevoflurane']['concentrations'], color='black', linestyle='-')
axes[3].set_title('Sevoflurane',fontsize=9)
axes[3].set_xlabel('Year',fontsize=9)
axes[3].set_ylabel('Annual Average (ppt)',fontsize=9)
axes[3].tick_params(axis='both', labelsize=7) 

# Adjust layout
plt.tight_layout()
plt.savefig('FigureS3')
plt.show()