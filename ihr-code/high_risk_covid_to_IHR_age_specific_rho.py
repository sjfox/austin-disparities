# -*- coding: utf-8 -*-
"""
This calculates the COVID IHR per age group and zip code, by first calculating
IHR for low and high risk individuals per age group based on national IHR for
France and the US, and then using high risk population numbers in each US zip
code.
"""
# %% Imports
import os
import pandas as pd
import numpy as np
# import itertools
# import copy
# import shapefile
#import matplotlib.pyplot as plt

np.set_printoptions(linewidth=125)
pd.set_option('display.width', 125)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# %% Functions
###############################################################################
def filter_df(df, conditions):
    # fltr = np.array([True]*len(df),dtype=bool)
    cond = conditions[0]
    col = cond[0]
    operation = cond[1]
    elts = cond[2]
    if operation == '==':
        fltr = (df[col].isin(elts))
    elif operation == '!=':
        fltr = ~(df[col].isin(elts))
    else:
        print('Wrong filtering operation in function filter_df')

    if len(conditions) > 1:
        for cond in conditions[1:]:
            col = cond[0]
            operation = cond[1]
            elts = cond[2]

            if operation == '==':
                fltr = fltr & (df[col].isin(elts))
            elif operation == '!=':
                fltr = fltr & ~(df[col].isin(elts))
            else:
                print('Wrong filtering operation in function filter_df')

    new_df = df.loc[fltr]

    return new_df
###############################################################################
def IHR_high_risk_calc(IHR_avg,prop_hr,rho):
    """ IHR_avg: IHR for the overall group
        prop_hr: proportion of high-risk individuals in that group. Number
    between 0 and 1
        rho: IHR(high-risk) / IHR(low-risk). Should be greater than 1.
        Calculate IHR for high-risk individuals as a function of the group's
    average IHR, the proportion of high-risk individuals, and the relative 
    IHR of high-risk to low-risk.
    """
    return IHR_avg / (prop_hr + (1 - prop_hr) / rho)
###############################################################################
# %% Load data
### Parameter: relative risk of high to low

### Input data
# Filenames and file locations
# folder_list = ['C:','Users','remyp','Research','COVID-19','High Risk']
# os.chdir(os.sep.join(folder_list))
folder_list = os.getcwd().split(os.sep)

population_folder = os.sep.join(folder_list + ['Population data','ZCTA',''])
population_filename = 'ZCTA Population per age group.xlsx'
pop_zcta_raw = pd.read_excel(population_folder + population_filename,sheet_name='Sheet1')

high_risk_zcta_filename = 'COVID High risk population per age group per zip code.csv'
high_risk_zcta_df = pd.read_csv(high_risk_zcta_filename)

rho_per_age_group_folder = os.sep.join(folder_list + ['Relative IHR per age-risk group',''])
rho_per_age_group_filename = 'Relative_IHR_age_risk_group.csv'
rho_per_age_group_df = pd.read_csv(rho_per_age_group_folder + rho_per_age_group_filename)
# rho = IHR(high-risk) / IHR(low-risk)

base_ihr = pd.read_excel('estimated-tx-ihr.xlsx',sheet_name='IHR')

texas_high_risk_filename = 'COVID High risk population per age group in Texas.csv'
texas_high_risk = pd.read_csv(texas_high_risk_filename)
df_fhr = texas_high_risk.iloc[:,1:]


# %% Keep ZCTA5 code only
pop_zcta = pop_zcta_raw.copy()
pop_zcta['ZCTA5'] = pop_zcta.apply(lambda x:
    x['GeographicAreaName'].replace('ZCTA5 ',''),axis=1)

## Population: split 0_4 into 0_0.5-4, and aggregate 75+
pop_zcta['0_0.5'] = pop_zcta['0_4']/8
pop_zcta['0.5_4'] = pop_zcta['0_4']*7/8

pop_zcta['75+'] = pop_zcta['75_79'] + pop_zcta['80_84'] + pop_zcta['85+']

del_cols = ['0_4','75_79','80_84','85+','GeographicAreaName','id']
for col in del_cols:
    if col in pop_zcta.columns:
        del pop_zcta[col]
        

## Rename population column: pop_ + age group
age_groups = [x for x in high_risk_zcta_df.columns if x != 'ZCTA5']
new_col_names = ['pop_' + x for x in age_groups]
pop_zcta.rename(columns=dict(zip(age_groups,new_col_names)),inplace=True)


# Format dataframes
high_risk_zcta_df['ZCTA5'] = high_risk_zcta_df['ZCTA5'].apply(str)
pop_zcta['ZCTA5'] = pop_zcta['ZCTA5'].apply(str)


###############################################################################
# %% Calculate IHR per zip code per age group, and average
## Get high-risk per age group in France or US in a dictionary

# Get rho per age group in a dictionary
rho_per_age_group_dict = dict(zip(
    rho_per_age_group_df['AgeGroup'], rho_per_age_group_df['rho']))
    
base_high_risk_d = dict(zip(df_fhr.columns,df_fhr.iloc[0]))

# Get IHR mean and bounds for high-risk individuals per age group
ihr_d = {}
ihr_high_risk_d = {}
for metric in ['Mean','LowerBound','UpperBound']:
    ihr_d[metric] = dict(zip(base_ihr['AgeGroup'],base_ihr[metric]))
    ihr_high_risk_d[metric] = {}
    for ag in age_groups:
        ihr_high_risk_d[metric][ag] = IHR_high_risk_calc(\
            ihr_d[metric][ag], base_high_risk_d[ag], rho_per_age_group_dict[ag])

# Scale US zip code numbers
# IHR_zcta = high_risk_zcta_df.copy()
IHR_zcta = None

for metric in ['Mean','LowerBound','UpperBound']:    # metric = 'Mean'
    IHR_m = high_risk_zcta_df.copy()
    IHR_m['Measure'] = metric
    
    for ag in age_groups:   #  ag = '20_24'
        ihr_hr = ihr_high_risk_d[metric][ag]
        IHR_m[ag] = ihr_hr * (IHR_m[ag] + (1-IHR_m[ag]) / rho_per_age_group_dict[ag])
        # IHR_m[ag] = np.min(IHR_m[ag],100)
    
    # Save results in the same dataframe
    if IHR_zcta is not None:
        IHR_zcta = pd.concat([IHR_zcta, IHR_m])
    else:
        IHR_zcta = IHR_m.copy()
         
        
## Calculate population weighted average per zip code
IHR_zcta = pd.merge(IHR_zcta,pop_zcta,on='ZCTA5',how='left')
IHR_zcta.rename(columns={'Total':'Population'},inplace=True)

IHR_zcta['AverageIHR'] = IHR_zcta.apply(lambda x:
    sum([x[y] * x['pop_' + y] for y in age_groups]) / x['Population']
    if (pd.notnull(x['Population']) and x['Population'] > 0) else np.nan, axis=1)

# Check a single zip code
# zip_code = '78741'
# IHR_zcta.loc[IHR_zcta['ZCTA5'] == zip_code,:]

out_cols = ['ZCTA5','Measure','Population','AverageIHR'] + age_groups 
IHR_zcta = IHR_zcta[out_cols]

out_filename_IHR_zcta = 'zip-age-ihr_estimated.csv'
IHR_zcta.to_csv(out_filename_IHR_zcta,index=False)
