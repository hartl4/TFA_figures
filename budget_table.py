#GENERATE BUDGET TABLE
# import modules
import xarray as xr
import numpy as np
import pandas as pd

#files
files_chem = 
files_tend = 

#Define Molecular weights
MW_AIR = 28.97
MW_TFA = 114
MW_TFAC = 132.4
MW_TFF = 116
MW_CF3CHO = 98


'''
#TERMS
PRODUCTION
tff_hyd_av: annual average Production of TFA from CF3COF
tfac_hyd_av: annual average Production of TFA from CF3COCl
cf3cho_hyd_av: annual average Production of TFA from CF3COH

LOSS
tff_loss_av: annual average loss of CF3COF through hydrolysis
tfac_loss_av: annual average loss of CF3COCl through hydrolysis
cf3cho_loss_av: annual average loss CF3COH through hydrolysis
tfac_phot_av = annual average loss of CF3COCl through photolysis
cf3cho_phot_av = annual average loss of CF3COH through photolysis
cf3cho_oh_av = annual average loss of CF3COH through OH loss
tfa_oh_av = annual average OH loss of TFA
tfac_dep_av = annual average wet + dry deposition of CF3COCl
tff_dep_av = annual average wet + dry deposition of CF3COF
cf3cho_dep_av = annual average wet + dry deposition of CF3COH
tfa_dep_av = annual average wet + dry deposition of TFA

BURDEN
tfac_burden_av =  annual average CF3COCl burden
tff_burden_av = annual average CF3COF burden
cf3cho_burden_av = annual average CF3COH burden
tfa_burden_av = annual average TFA burden

LIFETIME
tfac_life = lifetime of CF3COCl (days)
tff_life = lifetime of CF3COF (days)
cf3cho_life = lifetime of CF3COH (days)
tfa_life = lifetime of TFA (days)
'''
########################################################
#CALCULATE PRODUCTION OF TFA FROM CF3COX & LOSS OF CF3COX BY HYDROLYSIS

tfac_hyd = 0
tff_hyd = 0
cf3cho_hyd = 0

tfac_loss = 0
tff_loss = 0
cf3cho_loss = 0

years = 0

for file in files_chem:
    ds = xr.open_dataset(file,engine='netcdf4',decode_times=False)

    years += 1

    tfachyd = ds['chem191']
    tffhyd = ds['chem192']
    cf3chohyd = ds['chem193'] 

    #tfachyd = ds['chem128'] #SIMPLE
    #tffhyd = ds['chem129'] #SIMPLE
    #cf3chohyd = ds['chem130'] #SIMPLE
    
    #calculate production of tfa
    tfachyd_units = np.sum(tfachyd)*((MW_TFA/MW_AIR)*3600*1E-6) #Gg
    tffhyd_units = np.sum(tffhyd)*((MW_TFA/MW_AIR)*3600*1E-6) #Gg
    cf3chohyd_units = np.sum(cf3chohyd)*((MW_TFA/MW_AIR)*3600*1E-6) #Gg

    tfac_hyd += tfachyd_units
    tff_hyd += tffhyd_units
    cf3cho_hyd += cf3chohyd_units

    tfac_hyd_av = tfac_hyd/years
    tff_hyd_av = tff_hyd/years
    cf3cho_hyd_av = cf3cho_hyd/years

    #calculate loss of cf3cox
    tfacloss_units = np.sum(tfachyd)*((MW_TFAC/MW_AIR)*3600*1E-6) #Gg
    tffloss_units = np.sum(tffhyd)*((MW_TFF/MW_AIR)*3600*1E-6) #Gg
    cf3choloss_units = np.sum(cf3chohyd)*((MW_CF3CHO/MW_AIR)*3600*1E-6) #Gg

    tfac_loss += tfacloss_units
    tff_loss += tffloss_units
    cf3cho_loss += cf3choloss_units


    tfac_loss_av = tfac_loss/years
    tff_loss_av = tff_loss/years
    cf3cho_loss_av = cf3cho_loss/years

########################################################
#CALCULATE PRODUCTION OF CF3COX FROM SOURCE GASES

my_dict2 = {
    'tff': {"MW": 116E-3, "prod": ['chem082', 'chem178', 'chem091', 'chem181', 'chem097', 'chem183', 'chem117' ], "loss": 'chem192',  "yield":[]},
    'tfac': {"MW": 132.4E-3, "prod": ['chem078', 'chem175', 'chem086', 'chem179'], "loss": 'chem191',  "yield":[]},
    'cf3cho': {"MW": 98E-3, "prod": ['chem102', 'chem107',  'chem187', 'chem185', 'chem189'], "loss": ['chem119', 'chem193'], "yield":[]}} #base,base-min/max,base-lowres
'''
my_dict2 = {
    'tff': {"MW": 116E-3, "prod": ['chem077', 'chem079', 'chem080', 'chem084' ], "loss": 'chem129',  "yield":[]},
    'tfac': {"MW": 132.4E-3, "prod": ['chem076', 'chem078'], "loss": 'chem128',  "yield":[]},
    'cf3cho': {"MW": 98E-3, "prod": ['chem081', 'chem082',  'chem083'], "loss": ['chem086', 'chem130'], "yield":[]}}  #simple
'''

total_tfa_production = {key: 0 for key in my_dict2}  # Store total production for each precursor
num_years = len(files_chem)  # Total number of years in the dataset

for file in files_chem:
    ds = xr.open_dataset(file, engine='netcdf4', decode_times=False)
    year = file.split('em_')[1].split('_lat')[0]
    
    for key2 in my_dict2:
        total_prod_sum = 0
        
        # Convert to lists if needed
        prods = my_dict2[key2]["prod"]
        if not isinstance(prods, list):
            prods = [prods]

        # Calculate total production for each precursor
        for prod in prods:
            total_prod_sum += (ds[prod].sum() * 3600 * my_dict2[key2]["MW"] / 28.97E-3).values
        
        total_tfa_production[key2] += total_prod_sum  # Accumulate total production over years

# Compute the average annual TFA production for each precursor
avg_annual_tfa_production = {key: total_tfa_production[key]*1E-6 / num_years for key in my_dict2}

# Print the results
for key, avg_prod in avg_annual_tfa_production.items():
    print(f"Average annual TFA production from {key}: {avg_prod:.2f} Gg")

########################################################
#CALCULATE PHOTOLYSIS AND OH LOSS

tfac_phot = 0
cf3cho_phot = 0
tfa_oh = 0
cf3cho_oh = 0

years = 0

for file in files_chem:
    ds = xr.open_dataset(file,engine='netcdf4',decode_times=False)
    tfacphot = ds['chem171'] +ds['chem172'] #base
    cf3chophot = ds['chem173'] +ds['chem174'] #base
    ohloss = ds['chem076'] #base
    cf3cho_ohloss = ds['chem118'] #base

    tfacphot = ds['chem124'] +ds['chem125'] #SIMPLE
    cf3chophot = ds['chem126'] +ds['chem127'] #SIMPLE
    ohloss = ds['chem085'] #SIMPLE
    cf3cho_ohloss = ds['chem086'] #SIMPLE

    tfacphot_units = np.sum(tfacphot)*((MW_TFAC/MW_AIR)*3600*1E-6) #Gg
    cf3chophot_units = np.sum(cf3chophot)*((MW_CF3CHO/MW_AIR)*3600*1E-6) #Gg
    oh_units = np.sum(ohloss)*((MW_TFA/MW_AIR)*3600*1E-6) #Gg
    cf3cho_oh_units = np.sum(cf3cho_ohloss)*((MW_CF3CHO/MW_AIR)*3600*1E-6) #Gg

    tfac_phot += tfacphot_units
    cf3cho_phot += cf3chophot_units
    tfa_oh += oh_units
    cf3cho_oh += cf3cho_oh_units

    years += 1

    tfac_phot_av = tfac_phot/years
    cf3cho_phot_av = cf3cho_phot/years
    tfa_oh_av = tfa_oh/years
    cf3cho_oh_av = cf3cho_oh/years

########################################################
#CALCULATE DEPOSITION TERMS & STRATOSPHERIC LOSS

tfac_dep = 0
tff_dep = 0
cf3cho_dep = 0
tfa_dep = 0
tfac_strat = 0
tff_strat = 0
cf3cho_strat = 0
tfa_strat = 0
years = 0
tfa_wetdep = 0
tfa_drydep = 0

for file_tend in files_tend:
    ds = xr.open_dataset(file_tend,engine='netcdf4',decode_times=False)
    tfac_wet = ds['tfac_wetcnv']
    tfac_dry = ds['tfac_drydep']
    tff_wet = ds['tff_wetcnv']
    tff_dry = ds['tff_drydep']
    cf3cho_dry = ds['cf3cho_drydep']
    cf3cho_wet = ds['cf3cho_wetcnv']
    tfa_wet = ds['tfa_wetcnv']
    tfa_dry = ds['tfa_drydep']

    tff_strato = ds['tff_other']
    tfac_strato = ds['tfac_other']
    cf3cho_strato = ds['cf3cho_other']
    tfa_strato = ds['tfa_other']
    
    tfac = np.sum(tfac_wet + tfac_dry)
    tff = np.sum(tff_wet + tff_dry)
    cf3cho = np.sum(cf3cho_wet + cf3cho_dry)
    tfa = np.sum(tfa_wet + tfa_dry)

    tff_other = np.sum(tff_strato)
    tfac_other = np.sum(tfac_strato)
    cf3cho_other = np.sum(cf3cho_strato)
    tfa_other = np.sum(tfa_strato)
    
    tfac_dep += tfac
    tff_dep += tff
    cf3cho_dep += cf3cho
    tfa_dep += tfa
    tfa_wetdep += np.sum(tfa_wet)
    tfa_drydep += np.sum(tfa_dry)

    tfac_strat += tfac_other
    tff_strat += tff_other
    cf3cho_strat += cf3cho_other
    tfa_strat += tfa_other

    years += 1

    tfac_dep_av = tfac_dep/years*1E-6
    tff_dep_av = tff_dep/years*1E-6
    cf3cho_dep_av = cf3cho_dep/years*1E-6
    tfa_dep_av = tfa_dep/years*1E-6
    tfa_wet_av = tfa_wetdep/years*1E-6
    tfa_dry_av = tfa_drydep/years*1E-6

    tfac_strat_av = tfac_strat/years*1E-6
    tff_strat_av = tff_strat/years*1E-6
    cf3cho_strat_av = cf3cho_strat/years*1E-6
    tfa_strat_av = tfa_strat/years*1E-6

    
########################################################
#CALCULATE BURDENS FROM MASS FILES

tff_b = 0
tfac_b = 0
cf3cho_b = 0
tfa_b = 0 
count = 0

#open the netcdf file
for file_tend in files_tend:
    ds = xr.open_dataset(file_tend,engine='netcdf4',decode_times=False)
###calculate burdens from tend files
    air = ds['air_mass']
    tff= ds['tff_mass']
    tfac = ds['tfac_mass']
    cf3cho = ds['cf3cho_mass']
    tfa = ds['tfa_mass'] 

#convert VMR to MMR
    tfac_an = tfac.mean(dim='time')
    tff_an = tff.mean(dim='time')
    cf3cho_an = cf3cho.mean(dim='time')
    tfa_an = tfa.mean(dim='time')

    sum_tfac = np.sum(tfac_an)*1E-6
    sum_tff = np.sum(tff_an)*1E-6
    sum_cf3cho = np.sum(cf3cho_an)*1E-6
    sum_tfa = np.sum(tfa_an)*1E-6

    tfac_b += sum_tfac
    tff_b += sum_tff
    cf3cho_b += sum_cf3cho
    tfa_b += sum_tfa
    count +=1

    tfac_burden_av = tfac_b/count
    tff_burden_av = tff_b/count
    cf3cho_burden_av = cf3cho_b/count
    tfa_burden_av = tfa_b/count

########################################################
#CALCULATE TOTAL LOSS
tfac_loss_tot = tfac_phot_av + tfac_dep_av + tfac_loss_av + tfac_strat_av
cf3cho_loss_tot = cf3cho_phot_av + cf3cho_dep_av + cf3cho_loss_av + cf3cho_oh_av + cf3cho_strat_av
tff_loss_tot = tff_dep_av + tff_loss_av + tff_strat_av
tfa_loss_tot = tfa_dep_av + tfa_oh_av + tfa_strat_av

#CALCULATE TOTAL PRODUCTION
tfa_prod = tff_hyd_av + tfac_hyd_av + cf3cho_hyd_av

#CALCULATE LIFETIME
tfac_life = tfac_burden_av/tfac_loss_tot*365
tff_life = tff_burden_av/tff_loss_tot*365
cf3cho_life = cf3cho_burden_av/cf3cho_loss_tot*365
tfa_life = tfa_burden_av/tfa_loss_tot*365
########################################################

#PRINT TO EXCEL TABLE

def export_data():
    """
    Export atmospheric chemistry data to Excel table
    Replace the placeholder values with your actual calculated variables
    """
    
    # Helper function to extract values from xarrays
    def extract_value(var):
        """Extract numerical value from xarray or return as-is if already a number"""
        try:
            if hasattr(var, 'values'):
                # It's an xarray, extract the value
                val = var.values
                if hasattr(val, 'item'):  # Single value
                    return val.item()
                else:
                    return float(val)  # Convert to float
            else:
                return var  # Already a number
        except:
            return str(var)  # Fallback to string representation
    
    # Create the data structure
    data = {
        'Category': [],
        'Variable': [],
        'Description': [],
        'CF3COF': [],
        'CF3COCl': [],
        'CF3COH': [],
        'TFA': []
    }
    
    # PRODUCTION section
    production_data = [
        ('tff_hyd_av', 'Production of TFA from CF3COF', extract_value(tff_hyd_av), '-', '-', 'SOURCE'),
        ('tfac_hyd_av', 'Production of TFA from CF3COCl', '-', extract_value(tfac_hyd_av), '-', 'SOURCE'),
        ('cf3cho_hyd_av', 'Production of TFA from CF3COH', '-', '-', extract_value(cf3cho_hyd_av), 'SOURCE')
    ]
    
    for var, desc, cf3cof_val, cf3cocl_val, cf3coh_val, tfa_val in production_data:
        data['Category'].append('PRODUCTION')
        data['Variable'].append(var)
        data['Description'].append(desc)
        data['CF3COF'].append(cf3cof_val)
        data['CF3COCl'].append(cf3cocl_val)
        data['CF3COH'].append(cf3coh_val)
        data['TFA'].append(tfa_val)
    
    # LOSS section
    loss_data = [
        ('tff_loss_av', 'Loss through hydrolysis', extract_value(tff_loss_av), '-', '-', '-'),
        ('tfac_loss_av', 'Loss through hydrolysis', '-', extract_value(tfac_loss_av), '-', '-'),
        ('cf3cho_loss_av', 'Loss through hydrolysis', '-', '-', extract_value(cf3cho_loss_av), '-'),
        ('tfac_phot_av', 'Loss through photolysis', '-', extract_value(tfac_phot_av), '-', '-'),
        ('cf3cho_phot_av', 'Loss through photolysis', '-', '-', extract_value(cf3cho_phot_av), '-'),
        ('cf3cho_oh_av', 'Loss through OH reaction', '-', '-', extract_value(cf3cho_oh_av), '-'),
        ('tfa_oh_av', 'OH loss', '-', '-', '-', extract_value(tfa_oh_av))
    ]
    
    for var, desc, cf3cof_val, cf3cocl_val, cf3coh_val, tfa_val in loss_data:
        data['Category'].append('LOSS')
        data['Variable'].append(var)
        data['Description'].append(desc)
        data['CF3COF'].append(cf3cof_val)
        data['CF3COCl'].append(cf3cocl_val)
        data['CF3COH'].append(cf3coh_val)
        data['TFA'].append(tfa_val)
    
    # DEPOSITION section
    deposition_data = [
        ('tff_dep_av', 'Wet + dry deposition', extract_value(tff_dep_av), '-', '-', '-'),
        ('tfac_dep_av', 'Wet + dry deposition', '-', extract_value(tfac_dep_av), '-', '-'),
        ('cf3cho_dep_av', 'Wet + dry deposition', '-', '-', extract_value(cf3cho_dep_av), '-'),
        ('tfa_dep_av', 'Wet + dry deposition', '-', '-', '-', extract_value(tfa_dep_av)),
        ('tfa_dry_av', 'Dry deposition', '-', '-', '-', extract_value(tfa_dry_av)),
        ('tfa_wet_av', 'Wet deposition', '-', '-', '-', extract_value(tfa_wet_av))
    ]

    for var, desc, cf3cof_val, cf3cocl_val, cf3coh_val, tfa_val in deposition_data:
        data['Category'].append('DEPOSITION')
        data['Variable'].append(var)
        data['Description'].append(desc)
        data['CF3COF'].append(cf3cof_val)
        data['CF3COCl'].append(cf3cocl_val)
        data['CF3COH'].append(cf3coh_val)
        data['TFA'].append(tfa_val)
    # STRATOSPHERIC LOSS
    stratospheric_data = [
        ('tff_strat_av', 'Stratospheric loss', extract_value(tff_strat_av), '-', '-', '-'),
        ('tfac_strat_av', 'Stratospheric loss', '-', extract_value(tfac_strat_av), '-', '-'),
        ('cf3cho_strat_av', 'Stratospheric loss', '-', '-', extract_value(cf3cho_strat_av), '-'),
        ('tfa_strat_av', 'Stratospheric loss', '-', '-', '-', extract_value(tfa_strat_av))
    ]
    
    for var, desc, cf3cof_val, cf3cocl_val, cf3coh_val, tfa_val in stratospheric_data:
        data['Category'].append('DEPOSITION')
        data['Variable'].append(var)
        data['Description'].append(desc)
        data['CF3COF'].append(cf3cof_val)
        data['CF3COCl'].append(cf3cocl_val)
        data['CF3COH'].append(cf3coh_val)
        data['TFA'].append(tfa_val)
    
    # BURDEN section
    burden_data = [
        ('tff_burden_av', 'Annual average burden', extract_value(tff_burden_av), '-', '-', '-'),
        ('tfac_burden_av', 'Annual average burden', '-', extract_value(tfac_burden_av), '-', '-'),
        ('cf3cho_burden_av', 'Annual average burden', '-', '-', extract_value(cf3cho_burden_av), '-'),
        ('tfa_burden_av', 'Annual average burden', '-', '-', '-', extract_value(tfa_burden_av))
    ]
    
    for var, desc, cf3cof_val, cf3cocl_val, cf3coh_val, tfa_val in burden_data:
        data['Category'].append('BURDEN')
        data['Variable'].append(var)
        data['Description'].append(desc)
        data['CF3COF'].append(cf3cof_val)
        data['CF3COCl'].append(cf3cocl_val)
        data['CF3COH'].append(cf3coh_val)
        data['TFA'].append(tfa_val)
    
    # LIFETIME section (in days)
    lifetime_data = [
        ('tff_life', 'Lifetime (days)', extract_value(tff_life), '-', '-', '-'),
        ('tfac_life', 'Lifetime (days)', '-', extract_value(tfac_life), '-', '-'),
        ('cf3cho_life', 'Lifetime (days)', '-', extract_value(cf3cho_life), '-', '-'),
        ('tfa_life', 'Lifetime (days)', '-', extract_value(tfa_life), '-', '-')
    ]
    
    for var, desc, cf3cof_val, cf3cocl_val, cf3coh_val, tfa_val in lifetime_data:
        data['Category'].append('LIFETIME')
        data['Variable'].append(var)
        data['Description'].append(desc)
        data['CF3COF'].append(cf3cof_val)
        data['CF3COCl'].append(cf3cocl_val)
        data['CF3COH'].append(cf3coh_val)
        data['TFA'].append(tfa_val)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Export to Excel with formatting
    with pd.ExcelWriter('budget_table.xlsx', engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Results', index=False)
        
        # Get the workbook and worksheet
        workbook = writer.book
        worksheet = writer.sheets['Results']
        
        # Auto-adjust column widths
        for column in worksheet.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = min(max_length + 2, 50)
            worksheet.column_dimensions[column_letter].width = adjusted_width
    
    print("Excel file 'budget_table.xlsx' has been created successfully!")
    return df

export_data()  # For Excel
