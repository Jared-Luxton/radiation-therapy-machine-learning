import sys
import os
import pandas as pd
import numpy as np
from pandas import ExcelWriter
from pandas import ExcelFile
import re
from ast import literal_eval
import more_itertools
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib import lines
import imgkit
import seaborn as sns

from statistics import mean 
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from xgboost import XGBRegressor
from sklearn.metrics import explained_variance_score, r2_score
from sklearn.metrics import median_absolute_error
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import cross_val_score
import scipy.cluster.hierarchy as hac
import matplotlib.gridspec as gridspec
from scipy.stats import zscore
from scipy.stats import ks_2samp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
import scikit_posthocs as sp
from statsmodels.stats.anova import AnovaRM
       
def generate_dictionary_from_TeloLength_data(patharg):
    """
    opens raw telomere length measurement excel files from imageJ analyses and
    extracts the mean value of individual telomere length measurements to make histograms;
    telos are stored as values to their sample timepoint keys, which themselves 
    are values to patient id# key

    i.e the data structure is:

    dict = {
    patient_IDnumber: {SW#A non irrad: [telos data],
                       SW#A irrad @ 4 Gy: [telos data],
                       SW#B: [telos data],
                       SW#C: [telos data],
                       },
    etc.
    }

    i.e:

    all_patients_dict = {
    '1' = {'SW1A non irrad': ['telomere data'],
           'SW1A irrad @ 4 Gy': ['telomere data'],
           'SW1B': ['telomere data'],
           'SW1C': ['telomere data'],
          },
    etc. for patients 1 - 16 (less #4 missing)
    }
    USAGE:
    pass the directory where the telomere length excel files (.xlsx) are located
    """
    all_patients_dict = {}
    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
            try:
                df = pd.read_excel(file)
            except:
                print('File not found..')
                return -1

            print(file.name, 'data extraction in progress..')   
            telo_data = extract_and_clean_telos(df, file.name)
            file = file.name.replace('.xlsx', '').rstrip()
            data_list = []
            file_chr = ''
            num, num2 = capture_patient_sample_ID(file)

            if file[num:num2] not in all_patients_dict.keys():
                all_patients_dict[file[num:num2]] = {file: []}

                if len(all_patients_dict[file[num:num2]][file]) == 0:
                    all_patients_dict[file[num:num2]][file] = data_list
                    data_list.append(telo_data)
                    data_list.sort()
                elif len(all_patients_dict[file[num:num2]][file]) == 1:
                    data_list.append(telo_data)
                    data_list.sort()

            elif file[num:num2] in all_patients_dict.keys():
                if file in all_patients_dict[file[num:num2]]:
                    all_patients_dict[file[num:num2]][file].append(telo_data)
                    all_patients_dict[file[num:num2]][file].sort()
                elif file not in all_patients_dict[file[num:num2]]:     
                    all_patients_dict[file[num:num2]][file] = data_list
                    all_patients_dict[file[num:num2]][file].append(telo_data)
                    all_patients_dict[file[num:num2]][file].sort() 
                    
    print('completed file collection')
    return all_patients_dict


def generate_dataframe_from_dict(all_patients_dict):
    data = []
    for i in range(1,17):
        if str(i) in all_patients_dict.keys():
            for sample in sorted(all_patients_dict[str(i)].keys()):
                telos = all_patients_dict[str(i)][sample][0]
#                 IT WORKS PEGGY <333
                if 'hTERT' in sample:
                    #average of all hTERT samples is 79.9762
                    #CF = correction factor
                    hTERT_avg = 79.9762
                    hTERT_CF1 = hTERT_avg / telos['Mean Individ Telos'].mean()
                
                elif 'BJ1' in sample:
                    #average of all BJ1 samples is 69.5515
                    #CF = correction factor
                    BJ1_avg = 69.5515
                    BJ1_CF2 = BJ1_avg / telos['Mean Individ Telos'].mean()
   
                    CF_mean = (hTERT_CF1 + BJ1_CF2) / 2
            
                elif 'non irrad' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_A_nonRAD_name = sample
                    SW_A_nonRAD = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos)
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '1 ' + 'non irrad', telos_samp.multiply(CF_mean), individ_cells])

                elif 'irrad @ 4 Gy' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_A_irrad4Gy_name = sample
                    SW_A_irrad4Gy = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos)
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '2 ' + 'irrad @ 4 Gy', telos_samp.multiply(CF_mean), individ_cells])

                elif 'B' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_B_name = sample
                    SW_B = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos)
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '3 ' + 'B', telos_samp.multiply(CF_mean), individ_cells])
                    
                elif 'C' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_C_name = sample
                    SW_C = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos)
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '4 ' + 'C', telos_samp.multiply(CF_mean), individ_cells])

                else:
                    print('error with making dataframe from dict..')
                    print(sample)
                    continue
                     
    all_patients_df = pd.DataFrame(data, columns=['patient id', 'timepoint', 'telo data', 'cell data'])
    all_patients_df['patient id'] = all_patients_df['patient id'].astype('int')
    all_patients_df = all_patients_df.sort_values(by=['patient id', 'timepoint'], ascending=True, axis=0).reset_index(drop=True)
    all_patients_df['telo means'] = all_patients_df['telo data'].apply(lambda row: np.mean(row))
    
    all_patients_df['Q1'] = 'telos nonRAD Q1 <0.25'
    all_patients_df['Q2-3'] = 'telos nonRAD Q2-3 >0.25 & <0.75'
    all_patients_df['Q4'] = 'telos nonRAD Q4 >0.75'
    
    # counting telomeres per quartile
    all_patients_df = calculate_apply_teloQuartiles_dataframe(all_patients_df.copy())

    for quart in ['Q1', 'Q2-3', 'Q4']:
        all_patients_df[quart] = all_patients_df[quart].astype('int64')

    return all_patients_df


def gen_missing_values_andimpute_or_randomsampledown(n_cells, telosPercell, df):

    if df.size > 4600:
        dfsampled = df.sample(4600)
        return dfsampled

    if df.size > 25 and df.size <= 2300:
        missing_data_difference = abs((n_cells * telosPercell) - df.size)
        rsampled = df.sample(missing_data_difference, replace=True, random_state=28)
        rsampled = rsampled * 0.99999
        concat_ed = pd.concat([rsampled, df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed

    if df.size > 25 and df.size < 4600:
        missing_data_difference = abs((n_cells * telosPercell) - df.size)
        rsampled = df.sample(missing_data_difference, random_state=28)
        rsampled = rsampled * 0.99999
        concat_ed = pd.concat([rsampled, df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed
    
    else:
        print('not count standardized.. error')
        return df
    
    
def chunk_individual_telos_to_cells(telos_samp, n_telos):
    """
    splits up series of individual telomeres into equal parts, ala "cells"
    i.e if you have 92 telomeres per cell & 50 cells have contributed
    to this series, then pass 92 for n_telos.
    will return 50 cells each containing 92 telos 
    """
    telos_list = list(telos_samp)
    chunked_cells = more_itertools.chunked(telos_list, n_telos)
    chunked_cells = list(chunked_cells)
    
    cell_means = []
    
    for cell in chunked_cells:
        cell_means.append(np.mean(cell)) 
    
    return pd.Series(cell_means)
         

def capture_patient_sample_ID(file):
    if len(file) == 14:
        #it's patient id w/ 1 sample ID digit
        num = 2
        num2 = num + 1
        return num, num2
    
    elif len(file) == 12:
        #it's BJ1 ctrl w/ 1 sample ID digit
        num = 10
        num2 = num+ 1
        return num, num2
    
    elif 'hTERT' in file and len(file) == 17:
        #it's BJ-hTERT w/ 1 sample digit
        num = 15
        num2 = num + 1
        return num, num2
    
    elif len(file) == 15:
        #it's patient id w/ 2 sample ID digits
        num = 2
        num2 = num + 2
        return num, num2
    
    elif len(file) == 13:
        #it's BJ1 ctrl w/ 2 sample ID digits
        num = 10
        num2 = num + 2
        return num, num2
    
    elif 'hTERT' in file and len(file) == 18:
        # it's BJ-hTERT w/ 2 sample digits
        num = 15
        num2 = num + 2
        return num, num2
    
    elif len(file) == 4:
        #it's 2nd/3rd patient timepoint w/ 1 sample digit
        num = 2
        num2 = num + 1
        return num, num2
    
    elif len(file) == 5:
        #it's 2nd/3rd patient timepoint w/ 2 sample digits
        num = 2
        num2 = num + 2
        return num, num2
    
    elif '4 Gy' in file and len(file) == 17:
        # irrad @ 4 Gy 1 sample ID digit
        num = 2
        num2 = num + 1
        return num, num2
    
    elif '4 Gy' in file and len(file) == 18:
        # irrad @ 4 Gy 2 sample ID digits
        num = 2
        num2 = num + 2
        return num, num2

    
def extract_and_clean_telos(df, file_name):

    df.rename(columns={'Unnamed: 3':'Mean Individ Telos'}, inplace=True)
    mean_values_of_individual_telomere_lengths = (df['Mean Individ Telos'])
    mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.drop(
        labels=[5, 192, 379, 566, 753, 940, 1127, 1314, 1501, 1688, 1875, 2062, 2249, 2436, 2623, 2810, 2997, 3184, 
                3371, 3558, 3745, 3932, 4119, 4306, 4493, 4680, 4867, 5054, 5241, 5428, 5615, 5802, 5989, 6176, 6363, 
                6550, 6737, 6924, 7111, 7298, 7485, 7672, 7859, 8046, 8233, 8420, 8607, 8794, 8981, 9168])

    mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.iloc[6:9350]
    meantelos_str_toNaN = pd.to_numeric(mean_values_of_individual_telomere_lengths, errors='coerce')
    mean_individual_telos_cleaned = meantelos_str_toNaN.dropna(axis=0, how='any')
    mean_individ_df = mean_individual_telos_cleaned.to_frame(name=None)
    mean_individ_df.reset_index(drop=True, inplace=True)
    
    if 'BJ1' not in file_name and 'hTERT' not in file_name:
        telo_data = mean_individ_df[(np.abs(stats.zscore(mean_individ_df)) < 3).all(axis=1)]
        return telo_data
    else:
        return mean_individ_df
    

def quartile_cts_rel_to_df1(df1, df2):
    """
    FIND QUARTILES OF NON IRRAD TIMEPOINT & MAKE BASELINE..
    find individual telomeres below the 0.25 percentile (a), between
    the 0.25 & 0.75 percentile (b), & above the 0.75 percentile (c)
    """
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    
    quartile_1 = df2[df2 <= df1.quantile(0.25)].count()
    quartile_2_3 = df2[(df2 > df1.quantile(0.25)) & (df2 < df1.quantile(0.75))].count()
    quartile_4 = df2[df2 >= df1.quantile(0.75)].count()
    
    return float(quartile_1.values), float(quartile_2_3.values), float(quartile_4.values)


def calculate_apply_teloQuartiles_dataframe(all_patients_df):
    """
    LOOP THROUGH DATAFRAME FOR EACH PATIENT, ESTABLISH BASELINE QUARTILES FOR INDIVIDUAL TELOMERES USING NON IRRAD 
    SAMPLE TIMEPOINT.. THEN DETERMINES FOR EACH TIMEPOINT (irrad 4 Gy, B, C) HOW MANY TELOMERES REMAIN IN THOSE 
    QUARTILES... FILLS OUT Q1, Q2-3, Q4 COLUMNS..
    """
    
    q1_row, q2_3_row, q4_row = 5, 6, 7

    for i, row in all_patients_df.iterrows():
        if 'non irrad' in row[1]:
            nonRAD = row[2]
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] =                  (quartile_cts_rel_to_df1(nonRAD, nonRAD))
            
        elif 'irrad @ 4 Gy' in row[1]:
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, row[2]))

        elif 'B' in row[1]:
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, row[2]))

        elif 'C' in row[1]:
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, row[2]))

        else:
            print('unknown label in row[1] of the all patients df.. please check patient timepoint names')
            
    return all_patients_df


def histogram_plot_groups(x=None, data=None, groupby=None, iterable=None, n_bins=60, znorm=False):
    
    group_df = data.groupby(groupby)
    
    if groupby == 'timepoint':
            item = None
            non_irrad = group_df.get_group('1 non irrad').dropna(axis=0)[x]
            irrad_4_Gy = group_df.get_group('2 irrad @ 4 Gy').dropna(axis=0)[x]
            three_B = group_df.get_group('3 B').dropna(axis=0)[x]
            four_C = group_df.get_group('4 C').dropna(axis=0)[x]
            
            graph_four_histograms(non_irrad, n_bins, non_irrad, irrad_4_Gy, three_B, four_C,
                                                '1 non irrad', '2 irrad @ 4 Gy', '3 B', '4 C', znorm)
    
    elif groupby == 'patient id':
        for item in iterable:
            plot_df = group_df.get_group(item)
            
            non_irrad = plot_df[plot_df['timepoint'] == '1 non irrad'][x]
            irrad_4_Gy = plot_df[plot_df['timepoint'] == '2 irrad @ 4 Gy'][x]
            three_B = plot_df[plot_df['timepoint'] == '3 B'][x]
            four_C = plot_df[plot_df['timepoint'] == '4 C'][x]
            
            graph_four_histograms(non_irrad, n_bins, non_irrad, irrad_4_Gy, three_B, four_C,
                                  f'patient #{item} 1 non rad', 
                                  f'patient #{item} 2 irrad @ 4 Gy', 
                                  f'patient #{item} 3 B', 
                                  f'patient #{item} 4 C',
                                  znorm)
    #         plt.savefig(f'../graphs/telomere length/individual telos patient#{item}.png', dpi=400)

def graph_four_histograms(quartile_ref, n_bins, 
                          df1, df2, df3, df4,
                          name1, name2, name3, name4, znorm):
    
    n_bins = n_bins
    fig, axs = plt.subplots(2,2, sharey=True, sharex=True, 
                            constrained_layout=True, 
                            figsize = (8, 6))
    sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})

    histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df1, quartile_ref, name1, 0, 0, znorm)
    histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df2, quartile_ref, name2, 0, 1, znorm)
    histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df3, quartile_ref, name3, 1, 0, znorm)
    histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, df4, quartile_ref, name4, 1, 1, znorm)
    
#     plt.savefig(f'../graphs/telomere length/individual telos patient#{item}.png', dpi=400)

def histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone, axsNUMtwo, znorm):
    x_range = (0, 400)
    if znorm == True:
        x_range = (-3.5, 6.5) 
    astroarray = astroDF.to_numpy()
    N, bins, patches = axs[axsNUMone,axsNUMtwo].hist(astroarray, bins=n_bins, range=x_range, edgecolor='black')

    for a in range(len(patches)):
        if bins[a] <= np.quantile(astroquartile, 0.25):
            patches[a].set_facecolor('#fdff38')
        elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
            patches[a].set_facecolor('#d0fefe')
        elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
            patches[a].set_facecolor('#d0fefe')
        elif bins[a] > np.quantile(astroquartile, 0.75): 
            patches[a].set_facecolor('#ffbacd')
            
    axs[axsNUMone,axsNUMtwo].set_title(f"{astroname}", fontsize=14,)
    axs[axsNUMone,axsNUMtwo].tick_params(labelsize=12)
                
    font_axes=14

    if axsNUMone == 0 and axsNUMtwo == 0:
        axs[axsNUMone,axsNUMtwo].set_ylabel("Individual Telomere Counts", fontsize=font_axes)

    if axsNUMone == 1 and axsNUMtwo == 0:
        axs[axsNUMone,axsNUMtwo].set_ylabel("Individual Telomere Counts", fontsize=font_axes)
        axs[axsNUMone,axsNUMtwo].set_xlabel("Bins of Individual Telomeres (RFI)", fontsize=font_axes)

    if axsNUMone == 1 and axsNUMtwo == 1:
        axs[axsNUMone,axsNUMtwo].set_xlabel("Bins of Individual Telomeres (RFI)", fontsize=font_axes)

    axs[axsNUMone,axsNUMtwo].xaxis.set_major_locator(plt.MaxNLocator(12))
            
        


def color_seaborn_histogram(data, ax, n_bins, bin_vals):
    
    """
    rewriting individual telomere coloring scheme for seaborn.. from my original implementation in pandas/matplotlib

    provides access to values of bin edges:
    bin_vals = np.histogram(test, bins)[1]
    access to objects for coloring:
    ax.patches 
    
    usage:
    test = exploded_telos_all_patients_df[exploded_telos_all_patients_df['timepoint'] == '1 non irrad']['telo data exploded']
    ax = sns.set(font_scale=1)
    ax = sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
    bins = 80
    fig = plt.figure(figsize=(8,4))
    ax = sns.distplot(test, hist=True, kde=False, bins=bins, hist_kws=dict(alpha=.9))
    bin_vals = np.histogram(data, bins)[1]   
    
    color_seaborn_histogram(test, ax, bins)
    """   
        
    for a in range(len(ax.patches)):
        if bin_vals[a] < np.quantile(data, 0.25):
            ax.patches[a].set_facecolor('#fdff38')
        elif np.quantile(data, 0.25) < bin_vals[a] and bin_vals[a] <= np.quantile(data, 0.50):
            ax.patches[a].set_facecolor('#d0fefe')
        elif np.quantile(data, 0.50) < bin_vals[a] and bin_vals[a] <= np.quantile(data, 0.75):
            ax.patches[a].set_facecolor('#d0fefe')
        elif bin_vals[a] > np.quantile(data, 0.75): 
            ax.patches[a].set_facecolor('#ffbacd')
            
        
def plot_histogram_colored_stylizer(data=None, ax=None, n_bins=45):
    sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
    fig = plt.figure(figsize=(8,4))
    sns.distplot(data, hist=True, kde=False, bins=n_bins, hist_kws=dict(alpha=.9), ax=ax)
    bin_vals = np.histogram(data, n_bins)[1]   
    color_seaborn_histogram(data=data, ax=ax, n_bins=n_bins, bin_vals=bin_vals)
    
    
def grab_telo_data(patient_id_iterator=None, timepoint=None, df=None):
    data = df[(df['patient id'] == patient_id_iterator) &
              (df['timepoint'] == timepoint)]
    telos = data['individual telomeres']
    return telos
            
        
def make_timepoint_col(row):
    if 'A' in row:
        row = '1 non irrad'
        return row
    elif 'B' in row:
        row = '3 B'
        return row
    elif 'C' in row:
        row = '4 C'
        return row
    else:
        pass
        
        
def make_patient_ID(row):
    if len(row) == 4:
        return row[2]
    elif len(row) == 5:
        return row[2:4]


def change_sample_ID(row):
    if 'SW9C-2D' in row:
        row = 'SW9C'
        return row
    else:
        return row


#############################################################################

# Chromosome Aberration Methods 

#############################################################################

def make_dataframe_chr_aberr_data(patharg):
    all_chr_aberr_df = pd.DataFrame()
    for file in os.scandir(patharg):
            if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
                print(file)
                try:
                    df = pd.read_excel(file, usecols=list(range(29)), index_col=0, header=0)
                except:
                    print('File not found..')
                    return -1
            
                one_non_irrad = df.iloc[0:90]
                two_irrad_4_Gy = df.iloc[150:240]
                three_B = df.iloc[300:390]
                four_C = df.iloc[450:540]
                all_chr_aberr_df = pd.concat([one_non_irrad, two_irrad_4_Gy, three_B, four_C, all_chr_aberr_df])   
                
    return all_chr_aberr_df


def adjust_inversions_clonality(row):
    """
    df = df.apply(adjust_inversions_clonality, axis=1)
    """
    if row['sample notes'] == 'NaN' or row['sample notes'] == 'nan':
        pass
    if 'inv' in row['sample notes']:
        sample_notes = row['sample notes']
        clonal_inv = re.findall('[0-9] inv', sample_notes)
        
        if len(clonal_inv) > 0:
            row['# inversions'] = row['# inversions'] - len(clonal_inv)
        if 'term' in row['sample notes']:
            clonal_term_inv = re.findall('term inv', sample_notes)
    
            if len(clonal_term_inv) > 0:
                row['# terminal inversions'] = row['# terminal inversions'] - len(clonal_term_inv)
                
    if 'trans' in row['sample notes']:
        sample_notes = row['sample notes']
        clonal_trans = re.findall('[0-9] inv', sample_notes)
        
        if len(clonal_trans) > 0:
            row['translocations reciprocal 1,2,3'] = row['translocations reciprocal 1,2,3'] - len(clonal_trans)

    return row


#############################################
# MACHINE LEARNING HELPER FUNCTIONS / CLASSES


def stratify_SS_make_train_test(X, y, test_size, random_state):
    split = StratifiedShuffleSplit(n_splits=2, test_size=test_size, random_state=random_state)

    for train_index, test_index in split.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
            
    return X_train, X_test, y_train, y_test
    
    
def get_dummies_timepoints(df, timepts):
    col_retain_df = df[df['timepoint'].isin(timepts)]
    col_dummies = pd.get_dummies(col_retain_df, drop_first=True, columns=['timepoint'])
    return col_dummies


def calc_telomere_length_post_relative_pre(row):
    if row['4 C telo means'] > row['telo means']:
        row['length relative to pre'] = 1
    elif row['4 C telo means'] < row['telo means']:
        row['length relative to pre'] = 0
    return row 


def combine_data(exploded_telos=None, all_patients_df=None, 
                 prediction_objective='4 C means from individual telos',
                 timepoints_keep=['1 non irrad', '2 irrad @ 4 Gy']):
    
    if prediction_objective == '4 C means from individual telos': 
        col_to_rename = 'telo means'
        col_to_keep = 'individual telomeres'
        target = '4 C telo means'
        
    elif prediction_objective == '4 C means from telo means':
        col_to_rename = 'telo means'
        col_to_keep = 'telo means'
        target = '4 C telo means'
        
    elif prediction_objective == '4 C # short telos from individual telos':
        col_to_rename = 'Q1'
        col_to_keep = 'individual telomeres'
        target = '4 C # short telos'
    
    # pulling out 4 C
    four_C = all_patients_df[all_patients_df['timepoint'] == '4 C'][['patient id', col_to_rename, 'timepoint']]
    four_C.rename(columns={col_to_rename: target}, inplace=True)

    if prediction_objective == '4 C means from individual telos':
        # merging individual telomere data w/ 4 C telo means on patient id
        telo_data = (exploded_telos[exploded_telos['timepoint'] != '4 C']
                 .merge(four_C[[target, 'patient id']], on=['patient id']))
    
    elif prediction_objective == '4 C means from telo means':
        telo_data = (all_patients_df[all_patients_df['timepoint'] != '4 C']
             .merge(four_C[[target, 'patient id']], on=['patient id']))
    
    elif prediction_objective == '4 C # short telos from individual telos':
        telo_data = (exploded_telos[exploded_telos['timepoint'] != '4 C']
                 .merge(four_C[[target, 'patient id']], on=['patient id']))

    
#     cols to retain
#     cols_to_drop = [col for col in telo_data.columns if col not in cols_keep]
#     telo_data.drop(cols_to_drop, axis=1, inplace=True)
    telo_data = telo_data[['patient id', 'timepoint', col_to_keep, target]].copy()
    
    # timepoints of interest
    telo_data = telo_data[telo_data['timepoint'].isin(timepoints_keep)].copy()
    telo_data.reset_index(drop=True, inplace=True)
    return telo_data


class clean_data(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        return self
    
    
    def transform(self, X, y=None):       
        # renaming cols
        cols = list(X.columns)
        for col in cols:
            if ('timepoint_3 B' in cols) and ('_3 B' in col or 'irrad' in col):
                X.rename(columns={col:col.replace(' ', '')}, inplace=True)
            elif ('timepoint_3 B' not in cols) and ('irrad' in col):
                X.rename(columns={col:'timepoint'}, inplace=True)

        # enforcing col types
        cols = list(X.columns)
        for col in cols:
            if 'patient id' in col or 'timepoint' in col:
                X[col] = X[col].astype('int64')
            else:
                X[col] = X[col].astype('float64')

        return X


def grid_search(data, target, estimator, param_grid, scoring, cv, n_iter):
    
    grid = RandomizedSearchCV(estimator=estimator, param_distributions=param_grid, 
                              n_iter=n_iter, cv=cv, iid=False)
    
    pd.options.mode.chained_assignment = None  # this is because the gridsearch throws a lot of pointless warnings
    tmp = data.copy()
    grid = grid.fit(tmp, target)
    pd.options.mode.chained_assignment = 'warn'
    result = pd.DataFrame(grid.cv_results_).sort_values(by='mean_test_score', 
                                                        ascending=False).reset_index()
    
    del result['params']
    times = [col for col in result.columns if col.endswith('_time')]
    params = [col for col in result.columns if col.startswith('param_')]
    
    result = result[params + ['mean_test_score', 'std_test_score'] + times]
    
    return result, grid.best_estimator_


def cv_score_fit_mae_test(train_set=None, test_set=None, target='4 C telo means',
                          model=None, cv=5, scoring='neg_mean_absolute_error'):
    
    row = []
    features = [col for col in train_set if col != target and col != 'patient id']
    
    X_train = train_set[features].copy()
    X_test = test_set[features].copy()
    
    y_train = train_set[target].copy()
    y_test = test_set[target].copy()
    
    # cv
    scores = -1 * cross_val_score(model, X_train, y_train, cv=5, scoring=scoring)
    print(f'MAE per CV fold: \n{scores} \n')
    print(f'MEAN of MAE all folds: {scores.mean()}')
    print(f'STD of MAE all folds: {scores.std()}\n')

    # fitting the model
    model.fit(X_train, y_train)

    # predict y_test from X_test - this is using the train/test split w/o shuff;ing
    predict_y_test = model.predict(X_test)
    print(f"MAE of predict_y_test & y_test: {mean_absolute_error(y_test, predict_y_test)}")
    print(f'R2 between predict_y_test & y_test: {r2_score(y_test, predict_y_test)}')
    
    row.append(['XGBoost', features, target, round(scores.mean(), 4),
                                             round(scores.std(), 4),
                                             round(mean_absolute_error(y_test, predict_y_test), 4), 
                                             round(r2_score(y_test, predict_y_test), 4)])
    
    return model, row 


def predict_target_4C_compare_actual(telo_data=None, test_set=None, 
                                     model=None, target='4 C telo means',
                                     clean_process_pipe=None, verbose=True):

    telo_data = telo_data.copy()
    test_set_copy = test_set.copy()
    test_set_cleaner = clean_process_pipe
    if 'telo' in target:
        test_set_cleaned = test_set_cleaner.set_params(cleaner__drop_patient_id=False).fit_transform(test_set_copy)
    else:
        test_set_cleaned = test_set_cleaner.fit_transform(test_set_copy)
        
    features = [col for col in test_set if col != target]
    y_predict_list = []
    y_true_list = []

    for patient in list(telo_data['patient id'].unique()):
        # calculate actual mean telomere length per patient w/ all individual telos
        patient_data = telo_data[telo_data['patient id'] == patient]
        actual_4C = patient_data[target].mean()
        
        # calculate predicted mean telomere length per patient using only test data
        test_patient_data = test_set_cleaned[test_set_cleaned['patient id'] == patient].copy()
        test_patient_data.drop(['patient id', target], axis=1, inplace=True)
        predict_4C = model.predict(test_patient_data)
        if verbose:
            print(f'patient {patient}: ACTUAL {target}: {actual_4C:.2f} --- PREDICTED {target}: {np.mean(predict_4C):.2f}')
        y_predict_list.append(np.mean(predict_4C))
        y_true_list.append(actual_4C)
        
    print(f'MAE predicted vs. actual {target}: {mean_absolute_error(y_true_list, y_predict_list)}')
    print(f'R2 predicted vs. actual {target}: {r2_score(y_true_list, y_predict_list)}')
    
    return y_predict_list, y_true_list


class make_features(BaseEstimator, TransformerMixin):
    def __init__(self, make_log_individ_telos=False, make_log_target=False):
        self.make_log_individ_telos = make_log_individ_telos
        self.make_log_target = make_log_target
        
        
    def fit(self, X, y=None):
        return self
    
    
    def create_log_individ_telos(self, X, y=None):
        X['individual telos'] = np.log1p(X['individual telos'])
        return X
    
    
    def create_log_target(self, X, y=None):
        X['4 C telo means'] = np.log1p(X['4 C telo means'])
        return X
        
        
    def transform(self, X, y=None):
        if self.make_log_individ_telos:
            X = self.create_log_individ_telos(X)
            
        if self.make_log_target:
            X = self.create_log_target(X)
        return X
    
    
class make_dummies(BaseEstimator, TransformerMixin):
    def __init__(self, drop_first=True, cols_to_dummify=['timepoint']):
        self.drop_first = drop_first
        self.cols_to_dummify = cols_to_dummify
        
    
    def fit(self, X, y=None):
        return self
    
    
    def transf_dummies(self, X, y=None):
        dummies = pd.get_dummies(X, drop_first=self.drop_first, columns=self.cols_to_dummify)
        return dummies
    
    
    def transform(self, X, y=None):
        X = self.transf_dummies(X)
        return X
    
    
class clean_data(BaseEstimator, TransformerMixin):
    def __init__(self, drop_patient_id=True):
        self.drop_patient_id = drop_patient_id
    
    
    def fit(self, X, y=None):
        return self
    
    
    def transform(self, X, y=None):       
        # renaming cols
        cols = list(X.columns)
        i=1
        for col in cols:
            if ('timepoint_3 B' in cols) and ('_3 B' in col or 'irrad' in col):
                X.rename(columns={col:col.replace(' ', '')}, inplace=True)
#             elif ('timepoint_3 B' not in cols) and ('irrad' in col):
#                 X.rename(columns={col:f'timepoint_{i}'}, inplace=True)
#                 i+=1
                
        # enforcing col types
        cols = list(X.columns)
        for col in cols:
            if 'patient id' in col or 'timepoint' in col:
                X[col] = X[col].astype('int64')
            else:
                X[col] = X[col].astype('float64')
                
        if self.drop_patient_id:
            X.drop(['patient id'], axis=1, inplace=True)
            
        X.reset_index(drop=True, inplace=True)
        return X
    
    
########################################################################################################

######################               MACHINE LEARNING FOR CHR ABERR ...                 ################

########################################################################################################


class general_chr_aberr_cleaner(BaseEstimator, TransformerMixin):
    def __init__(self, adjust_clonality=True, combine_alike_aberr=True, drop_what_timepoint='3 B'):
        self.adjust_clonality = adjust_clonality
        self.combine_alike_aberr = combine_alike_aberr
        self.drop_what_timepoint = drop_what_timepoint
        
    def fit(self, X, y=None):
        return self
    
    
    def enforce_column_types(self, X, y=None):
        for col in X.columns:
            if col == 'sample notes' or col == 'timepoint':
                X[col] = X[col].astype('str')
            elif col in ['# inversions', '# terminal inversions', 
                         'translocations reciprocal 1,2,3',
                         'translocations one-way 1,2,3', 'dicentrics']:
                X[col] = X[col].astype('int64')
        return X
    
    
    def adjust_clonality_counts(self, row):
        if row['sample notes'] == 'NaN' or row['sample notes'] == 'nan':
            pass
    
        if 'inv' in row['sample notes'] and 'term' not in row['sample notes']:
            sample_notes = row['sample notes']
            clonal_inv = re.findall('[0-9] inv', sample_notes)
#             print(sample_notes, clonal_inv)

            if len(clonal_inv) > 0:
                row['# inversions'] = row['# inversions'] - len(clonal_inv)
                
        if 'term' in row['sample notes']:
            sample_notes = row['sample notes']
            clonal_term_inv = re.findall('[0-9] term', sample_notes)
#             print(sample_notes, clonal_term_inv)
                
            if len(clonal_term_inv) > 0:
                row['# terminal inversions'] = row['# terminal inversions'] - len(clonal_term_inv)

        if 'trans' in row['sample notes']:
            sample_notes = row['sample notes']
            clonal_trans = re.findall('[0-9] trans', sample_notes)
            
            if len(clonal_trans) > 0:
                row['translocations reciprocal 1,2,3'] = row['translocations reciprocal 1,2,3'] - len(clonal_trans)

        return row

    
    def combine_chr_aberr(self, X, y=None):
        # combining satellite associations
        X['# sat associations'] = (X['# 2 chr sat. associations'] + X['# 3 chr sat. associations'] +
                                   X['# 4 chr sat. associations'] + X['# 5+ chr sat. associations'])
        # combining terminal SCEs
        X['# terminal SCEs'] = X['# terminal SCEs cis-paint'] + X['# terminal SCEs cis-dark']
        # combining translocations
        X['# translocations'] = X['translocations reciprocal 1,2,3'] + X['translocations one-way 1,2,3']
        return X
    
    
    def drop_columns(self, X, y=None):
        # dropping unneeded chr aberr types
        X = X.drop(columns=['# 2 chr sat. associations', '# 3 chr sat. associations', 
                            '# 4 chr sat. associations', '# 5+ chr sat. associations',
                            '# terminal SCEs cis-paint', '# terminal SCEs cis-dark',
                            'translocations reciprocal 1,2,3', 'translocations one-way 1,2,3',
                            'tricentrics', '# sub-telo SCEs',
                            'chr fragments', 'expected chr fragments'
                           ], axis=1)
        # dropping misc. notation columns
        X = X.drop(columns=['metaphase size', 'terminal inversion size', 'inversion size',
                            'inversion notes', 'terminal inversion notes',
                            'translocation intra notes', 
                            'sample notes', 
                            'chromosome'
                           ], axis=1)
        return X
    
    
    def drop_timepoint(self, X, y=None):
        X = X[X['timepoint'] != self.drop_what_timepoint].copy()
        return X
    
    
    def drop_patient_ID(self, X, y=None):
        X = X[X['patient id'] != 13].copy()
        return X

    
    def transform(self, X, y=None):
        X = self.enforce_column_types(X)
        
        if self.adjust_clonality:
            X = X.apply(self.adjust_clonality_counts, axis=1)
        if self.combine_alike_aberr:
            X = self.combine_chr_aberr(X)
        if self.drop_what_timepoint:
            X = self.drop_timepoint(X)
        
        X = self.drop_columns(X)
        X = self.drop_patient_ID(X)
        
        X.rename(columns={'dicentrics': '# dicentrics'}, inplace=True)
        # data is arranged as events (chr aberr) per chromosome per cell; first sum per cell
        X = X.groupby(['patient id', 'timepoint', 'cell number']).agg('sum').reset_index()
        X.drop(['cell number'], axis=1, inplace=True) 
        return X
    
    
class make_chr_features(BaseEstimator, TransformerMixin):
    def __init__(self, combine_inversions=False, bool_features=False,
                       features=['# inversions', '# terminal inversions', '# dicentrics', '# translocations']):
        self.combine_inversions = combine_inversions
        self.bool_features = bool_features
        self.features = features
    
    
    def fit(self, X, y=None):
        return self
    
    
    def total_inversions(self, X, y=None):
        X['# total inversions'] = X['# inversions'] + X['# terminal inversions']
        return X
    
    
    def true_false(self, row):
        if int(row) == 0:
            row = False
        elif int(row) > 0:
            row = True
        return row
    
    def make_bool_features(self, X, y=None):
        for feature in self.features:
            X[f'BOOL {feature}'] = X[feature].apply(lambda row: self.true_false(row))
            X[f'BOOL {feature}'] = X[f'BOOL {feature}'].astype('bool')
        return X
    
    
    def transform(self, X, y=None):
        if self.combine_inversions:
            X = self.total_inversions(X)
        if self.bool_features:
            X = self.make_bool_features(X)
        return X
    
    
class make_target_merge(BaseEstimator, TransformerMixin):
    def __init__(self, target='aberration index', target_timepoint='4 C', target_type='means',
                       features=['# inversions', '# terminal inversions', '# dicentrics', '# translocations'], 
                       drop_first=True):
        self.target = target
        self.target_timepoint = target_timepoint
        self.target_type = target_type
        self.features = features
        self.drop_first = drop_first

        
    def fit(self, X, y=None):
        return self
    
    
    def encode_target_4C(self, row):
        target = self.target
        irrad4gy_index = '2 irrad @ 4 Gy mean aberration index'
        encoded = '4 C encoded mean aberration index'
        
        if row[target] > row[irrad4gy_index]:
            row[encoded] = 2
        elif row[target] == row[irrad4gy_index]:
            row[encoded] = 1
        elif row[target] < row[irrad4gy_index]:
            row[encoded] = 0
        else:
            print('error')
        return row
    
    
    def extract_timepoint_rename(self, X, y=None, timept=None):
        target_col = f'{timept} {self.target}'
        timept_means = X[X['timepoint'] == timept][['patient id', 'timepoint', self.target]]
        timept_means.rename(columns={self.target: target_col}, inplace=True)
        return timept_means, target_col
    
    
    def arrange_features_target(self, X, y=None):
        bool_cols = [col for col in X.columns if 'BOOL' in col]
        X = X[['patient id', 'timepoint'] + self.features + bool_cols].copy()
        if self.target == 'aberration index':
            X[self.target] = X[self.features].sum(axis=1)
            
        X_means = X.groupby(['patient id', 'timepoint']).agg('mean').reset_index()
        if self.target == 'aberration index':
            X.drop([self.target], axis=1, inplace=True)
        fourC_means, target_4C = self.extract_timepoint_rename(X_means, timept='4 C')
#         irrad4Gy_means, target_irr4Gy = self.extract_timepoint_rename(X_means, timept='2 irrad @ 4 Gy')
        
        complete = X.merge(fourC_means[['patient id', target_4C]], on='patient id')
#         complete = complete.merge(irrad4Gy_means[['patient id', target_irr4Gy]], on='patient id')
        complete = complete[complete['timepoint'] != '4 C'].copy()        
#         complete.drop([target_irr4Gy], axis=1, inplace=True)
        return complete
        
    
    def encode_timepoint_col(self, X, y=None):
        dummies = pd.get_dummies(X, drop_first=self.drop_first, columns=['timepoint'])
        return dummies
    
     
    def transform(self, X, y=None):
        X = self.arrange_features_target(X)
        X = self.encode_timepoint_col(X)
        return X
    
    
def xgboost_hyper_param(learning_rate, n_estimators, max_depth,
                        subsample, colsample, gamma):
 
    max_depth = int(max_depth)
    n_estimators = int(n_estimators)
 
    clf = XGBRegressor(max_depth=max_depth,
                       learning_rate=learning_rate,
                       n_estimators=n_estimators,
                       gamma=gamma, objective='reg:squarederror')
    
    return np.mean(cross_val_score(clf, X_train, y_train, cv=5, scoring='neg_mean_absolute_error'))


    encode_dict = {'1 non irrad' : 1, '2 irrad @ 4 Gy': 2, '3 B': 3, '4 C': 4}
    return encode_dict[row]


def chr_aberr_predict_target_4C_compare_actual(cleaned_unsplit_chr_data=None, cleaned_test_set=None, 
                                     model=None, target='4 C aberration index',
                                     clean_process_pipe=None, verbose=True):

    chr_data = cleaned_unsplit_chr_data.copy()
    chr_data = clean_process_pipe.fit_transform(chr_data)
    chr_test = cleaned_test_set.copy()
#     test_set_cleaned = clean_process_pipe.fit_transform(test_set_copy)
        
    features = [col for col in chr_test if col != target]
    y_predict_list = []
    y_true_list = []

    for patient in list(chr_data['patient id'].unique()):
        # calculate actual mean telomere length per patient w/ all individual telos
        patient_data = chr_data[chr_data['patient id'] == patient]
#         display(patient_)
        actual_4C = patient_data[target].mean()
        
        # calculate predicted mean telomere length per patient using only test data
        test_patient_data = chr_test[chr_test['patient id'] == patient].copy()
        test_patient_data.drop(['patient id', target], axis=1, inplace=True)
        predict_4C = model.predict(test_patient_data)
        if verbose:
            print(f'patient {patient}: ACTUAL {target}: {actual_4C:.2f} --- PREDICTED {target}: {np.mean(predict_4C):.2f}')
        y_predict_list.append(np.mean(predict_4C))
        y_true_list.append(actual_4C)
        
    print(f'MAE predicted vs. actual {target}: {mean_absolute_error(y_true_list, y_predict_list)}')
    print(f'R2 predicted vs. actual {target}: {r2_score(y_true_list, y_predict_list)}')
    
    return y_predict_list, y_true_list


################################################################################

###############              CLUSTERING ANALYSES                 ###############

################################################################################


def encode_timepts(row):
    encode_dict = {'1 non irrad' : 1, '2 irrad @ 4 Gy': 2, '3 B': 3, '4 C': 4}
    return encode_dict[row]


def myMetric(x, y):
    r = stats.pearsonr(x, y)[0]
    return 1 - r 


def plot_dendogram(Z, target=None, indexer=None):
    with plt.style.context('fivethirtyeight' ): 
        plt.figure(figsize=(10, 2.5))
        plt.title(f'Dendrogram of clusters by {target}', fontsize=22, fontweight='bold')
        plt.xlabel('patient IDs', fontsize=22, fontweight='bold')
        plt.ylabel('distance', fontsize=22, fontweight='bold')
        hac.dendrogram(Z, labels=indexer, leaf_rotation=90.,    # rotates the x axis labels
                        leaf_font_size=15., ) # font size for the x axis labels
        plt.show()

        
def plot_results(timeSeries, D, cut_off_level, y_size, x_size, verbose):
    result = pd.Series(hac.fcluster(D, cut_off_level, criterion='maxclust'))
    if verbose:
        clusters = result.unique() 
        fig = plt.subplots(figsize=(x_size, y_size))   
        mimg = math.ceil(cut_off_level/2.0)
        gs = gridspec.GridSpec(mimg,2, width_ratios=[1,1])   
        for ipic, c in enumerate(clusters):
            cluster_index = result[result==c].index
            print(ipic, "Cluster number %d has %d elements" % (c, len(cluster_index)))
            ax1 = plt.subplot(gs[ipic])
            ax1.plot(timeSeries.T.iloc[:,cluster_index])
            ax1.set_title((f'Cluster number {c}'), fontsize=15, fontweight='bold')      
        plt.show()
    return result
        
        
def cluster_data_return_df(df, target='telo means', cut_off_n=4, 
                           metric=myMetric, method='single',
                           y_size=6, x_size=10, verbose=True):
    
    df = df[df['patient id'] != 13].copy()
    # preparing data
    if '1 non irrad' in df['timepoint'].unique():
        df['timepoint'] = df['timepoint'].apply(lambda row: encode_timepts(row))
    df = df[['patient id', 'timepoint', target]].copy()
    df = df.pivot(index='patient id', values=target, columns='timepoint').reset_index()
#     if 13 in df['patient id'].unique() and 'telo means' in df.columns:
#         df.drop(11, inplace=True, axis=0)
    df.set_index('patient id', inplace=True)
    
    # run the clustering    
    cluster_Z = hac.linkage(df, method=method, metric=metric)
    if verbose:
        plot_dendogram(cluster_Z, target=target, indexer=df.index)
    # return df bearing cluster groups
    indexed_clusters = plot_results(df, cluster_Z, cut_off_n, y_size=y_size, x_size=x_size, verbose=verbose)
    
    # concat clusters to original df and return
    ready_concat = df.reset_index()
    clustered_index_df = pd.concat([ready_concat, indexed_clusters], axis=1)
    clustered_index_df.rename(columns={clustered_index_df.columns[-1]: f'{target} cluster groups',
                                       1: '1 non irrad',
                                       2: '2 irrad @ 4 Gy',
                                       3: '3 B',
                                       4: '4 C'}, inplace=True)
    melted = clustered_index_df.melt(id_vars=['patient id', f'{target} cluster groups'], 
                                     var_name='timepoint', value_name=target)
    return melted


def graph_cluster_groups(df, target=None, hue=None):
    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    
    ax = sns.set(font_scale=1.5)
    ax = sns.lineplot(x='timepoint', y=target, data=df, hue=hue, legend='full',
            palette=sns.color_palette(flatui[:len(df[hue].unique())]))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.setp(ax.get_xticklabels(), rotation=45)
    

def graph_clusters_per_patient(df, target=None, cluster_name=None,
                               y_dimen=2, x_dimen=2, fsize=(9,8)):
    if cluster_name == None:
        cluster_name = f'{target} cluster groups'
        
    fig, ax = plt.subplots(y_dimen,x_dimen, sharex='col', sharey='row', figsize=fsize)
    axes = ax.ravel()
    n_groups = df[cluster_name].nunique()

    for i in range(1, n_groups + 1):
        data_clusters = df[df[cluster_name] == i]
        print(f'{target} CLUSTER {i} | patient IDs: {list(data_clusters["patient id"].unique())}')
        sns.lineplot(x='timepoint', y=target, data=data_clusters, hue='patient id', legend=False, ax=axes[i-1],
                     palette=sns.color_palette("Set1", data_clusters['patient id'].nunique()))
        axes[i-1].set_title((f'Cluster number {i}'), fontsize=15, fontweight='bold')      
    for ax in fig.axes:
        plt.setp(ax.get_xticklabels(), horizontalalignment='right', rotation=45)
              
              
def eval_make_test_comparisons(df=None, timepoints=None, test=None, test_name=None, 
                               target='individual telomeres'):
    if timepoints == None:
        timepoints = list(df['timepoint'].unique())
    timept_pairs = []
    row = []

    g_1 = df[df['timepoint'] == '1 non irrad'][target]
    g_2 = df[df['timepoint'] == '2 irrad @ 4 Gy'][target]
    g_3 = df[df['timepoint'] == '3 B'][target]
    g_4 = df[df['timepoint'] == '4 C'][target]
    df_list = [g_1, g_2, g_3, g_4]
              
    for iter1, df in zip(timepoints, df_list):
        for iter2, i in zip(timepoints, range(len(df_list))):
            pair1, pair2 = f"{iter1}:{iter2}", f"{iter2}:{iter1}"
            if iter2 != iter1 and pair1 not in timept_pairs and pair2 not in timept_pairs:
                stat, pvalue = test(df, df_list[i])
                print(f'{test_name} | {iter1} vs {iter2} {pvalue}')
                timept_pairs.append(pair1)
                timept_pairs.append(pair2)
                row.append([test_name, iter1, iter2, pvalue])
    return timept_pairs, row
              
              
################################################################
################                                ################
################             MISC               ################
################                                ################            
################################################################

              
def telos_scipy_anova_post_hoc_tests(df=None, time_col='timepoint', target='individual telomeres',
                                     sig_test=stats.f_oneway, post_hoc=None):

    g_1 = df[df['timepoint'] == '1 non irrad'][target]
    g_2 = df[df['timepoint'] == '2 irrad @ 4 Gy'][target]
    g_3 = df[df['timepoint'] == '3 B'][target]
    g_4 = df[df['timepoint'] == '4 C'][target]
    statistic, p_value = sig_test(g_1, g_2, g_3, g_4)
    print(f'ONE WAY ANOVA for telomere length: {p_value}')
              
    # if anova detects sig diff, perform post-hoc tests
    if p_value <= 0.05:
        mc = MultiComparison(df[target], df['timepoint'])
        mc_results = mc.tukeyhsd()
        print(mc_results)
              
              
def chr_scipy_anova_post_hoc_tests(df=None, flight_status_col='timepoint',
                                   sig_test=stats.f_oneway, post_hoc=sp.posthoc_ttest):
    """
    df should be melted by aberration type
    """
    # make list of aberrations
    aberrations = list(df['aberration type'].unique())
    
    # loop through aberrations & perform anovas between pre/mid/post
    for aberr in aberrations:
        g_1 = df[(df[flight_status_col] == '1 non irrad') & (df['aberration type'] == aberr)]['count per cell']
        g_2 = df[(df[flight_status_col] == '2 irrad @ 4 Gy') & (df['aberration type'] == aberr)]['count per cell']
        g_3 = df[(df[flight_status_col] == '3 B') & (df['aberration type'] == aberr)]['count per cell']
        g_4 = df[(df[flight_status_col] == '4 C') & (df['aberration type'] == aberr)]['count per cell']
        statistic, p_value = sig_test(g_1, g_2, g_3, g_4)
        print(aberr, p_value)
        
        # if anova detects sig diff, perform post-hoc tests
        if p_value <= 0.05:
            if post_hoc == 'tukeyHSD':
                mc = MultiComparison(df[df['aberration type'] == aberr]['count per cell'], 
                                     df[df['aberration type'] == aberr]['timepoint'])
                print(mc.tukeyhsd())
            else:
                display(post_hoc(df[df['aberration type'] == aberr], val_col='count per cell', 
                        group_col='timepoint', equal_var=False))
        print('\n')
              

def z_norm_individual_telos(exploded_telos_df=None):
    z_norm = pd.DataFrame()
    df = exploded_telos_df
    grouped_df = df.groupby('timepoint')
    
    for timept in df['timepoint'].unique():
        timept_df = grouped_df.get_group(timept).copy()
        timept_df['z-norm_individual_telos'] = zscore(timept_df['individual telomeres'])
        timept_df.reset_index(drop=True, inplace=True)
        z_norm = pd.concat([z_norm, timept_df], axis=0)
    return z_norm
              
              
def script_load_clean_data_ml_pipeline_loop_aberrations(features_list=None, target1_list=None, target2_list=None, stats_list=None,
                                                        verbose=True):
    graphing_dict = {}
    for features, target1, target2 in zip(features_list, target1_list, target2_list):
        # loading chr aberr data
        all_chr_aberr_df = pd.read_csv('../compiled patient data csv files/all_chr_aberr_df.csv')
        # initializing general cleaner pipeline for data
        general_cleaner = Pipeline([('cleaner', general_chr_aberr_cleaner())])
        # cleaning data
        cleaned_chr_df = general_cleaner.fit_transform(all_chr_aberr_df)
        #splitting cleaned data into train/test
        chr_train, chr_test = train_test_split(cleaned_chr_df, test_size=0.2, shuffle=True, 
                                           stratify=cleaned_chr_df[['patient id', 'timepoint']])
        # initializing pipeline to generate features + target from data for machine learning
        make_new_features_target = Pipeline([('make features', make_chr_features(combine_inversions=False, bool_features=False, 
                                                                                     features=features)),
                                             ('make target merge', make_target_merge(target=target1, features=features))])
        # making new train/test dataframes w/ features & target
        cleaned_chr_train = chr_train.copy()
        cleaned_chr_test = chr_test.copy()
        cleaned_chr_train = make_new_features_target.fit_transform(cleaned_chr_train)
        cleaned_chr_test = make_new_features_target.fit_transform(cleaned_chr_test)

        # creating XGBoost regressor model 
        chr_model = XGBRegressor(n_estimators=200, max_depth=15, learning_rate=0.1, objective='reg:squarederror', random_state=0,)

        # performing cross-fold validation of XGBoost regressor on training set, returns model fitted on training data
        print(f'--------------------------------------------------------------------')
        print(f'PERFORMING CROSSFOLD VALIDATION of XGBoost model') 
        print(f'features: {features} ___ target: {target2}')
        print(f'--------------------------------------------------------------------')
        chr_fit_xgb_model, row = cv_score_fit_mae_test(train_set=cleaned_chr_train, test_set=cleaned_chr_test,
                                                           model=chr_model, cv=5, target=target2)
        stats_list += row

        # predicting target from test data w/ trained model & comparing predicted vs. actual values
        print(f'--------------------------------------------------------------------')
        print(f'PREDICTIONS of trained XGBoost model vs. actual values') 
        print(f'features: {features} ___ target: {target2}')
        print(f'--------------------------------------------------------------------')
        chr_y_predict, y_true = chr_aberr_predict_target_4C_compare_actual(cleaned_unsplit_chr_data=cleaned_chr_df,
                                                                           cleaned_test_set=cleaned_chr_test, 
                                                                           model=chr_fit_xgb_model, target=target2,
                                                                           clean_process_pipe=make_new_features_target,
                                                                           verbose=verbose)
        graphing_dict[target1] = [y_true, chr_y_predict]
        print('\n')
        
    return stats_list, graphing_dict
              
              
def make_stats_df(stats_list=None):    
    stats_df = pd.DataFrame(data=stats_list, columns=['Model', 'Features', 'Target', 
                                                      'Average MAE of CV folds', 'Std dev of MAE of CV folds', 
                                                      'MAE predicted vs. test values', 'R2 predicted vs. test values'])
    return stats_df
              
def make_graphing_df(graphing_dict=None):
    graphing_df = pd.DataFrame()
    for key in graphing_dict.keys():
        data = pd.DataFrame({'aberration type':key, 
                             'actual values':graphing_dict[key][0], 
                             'predicted values':graphing_dict[key][1]})
        graphing_df = pd.concat([graphing_df, data], axis=0)
    return graphing_df
              
def plot_individ_telos_ML_objective(df=None, timept_col='timepoint', 
                                    timept_1='1 non irrad', timept_2='2 irrad @ 4 Gy',
                                    features='individual telomeres',
                                    target='4 C telo means'):
    # create subplot object
    fig, ax = plt.subplots(3, 5, figsize=(24,15), 
#                            sharex='col', 
                           sharey='row'
                          )

    # create flattened axes, loop through axes and populate w/ histograms & data
    axes = ax.ravel()
    for i, id_num in zip(axes, df['patient id'].unique()):
        df2 = df[df['patient id'] == id_num].copy()
        nonirrad = df2[df2[timept_col] == timept_1][features]
        irrad4Gy = df2[df2[timept_col] == timept_2][features]

        # plotting 1 non irrad & 2 irrad @ 4 Gy distributions + ML target
        sns.distplot(nonirrad, bins=25, hist=True, ax=i, label=timept_1, norm_hist=False, kde=False)
        sns.distplot(irrad4Gy, bins=25, hist=True, ax=i, label=timept_2, norm_hist=False, kde=False)
        i.axvline(x=df2[target].unique(), c='black', linewidth=3, alpha=0.5, linestyle='--')
        
        # labeling subplots
        i.tick_params(labelsize=14)
        i.set_xlabel('')
        i.set_title(f'patient #{id_num}', fontsize=18)
    axes[-1].axis('off')
#     axes[-1].grid(False)
#     axes[-1].set_facecolor('white')
#     axes[-1].set_yticks([])
#     axes[-1].set_title('')

    # create legend
    handles, labels = i.get_legend_handles_labels()
    # create line object to represent 4 C target drawn on plots & append to legend objects
    vertical_line = lines.Line2D([], [], color='black', linestyle='--')
    handles.append(vertical_line)
    labels.append(f'{target} per patient')
    plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=22)

    # label axes
    fig.text(0.5, -.005, 'Bins of Individual Telomeres (RFI)', ha='center', va='center', fontsize=24)
    fig.text(-.005, 0.5, 'Counts of Individual Telomeres', ha='center', va='center', rotation='vertical', fontsize=24)
    plt.tight_layout()
    
    # save
    plt.savefig(f'../graphs/paper figures/supp figs/visualize {target} objective individual telos ML model.png', 
                dpi=400, bbox_inches = "tight")
              
              
def df_to_png(df=None, path=None):
    html_df = df.to_html(justify='center')
    imgkit.from_string(html_df, path, options = {'format': 'png'})
              
                      
def plot_multiple_types_clusters(y_list=None, hue_list=None,
                                 df_list=None, ylim_dict=None):
              
    # defining plot size, layout based on num clusters to plot
    fsize=(9.6,4.8)
    row_dim = 1
    col_dim = 2     
    scale=1.5
    if len(y_list) > 2:
        fsize=(13.6,8)
        row_dim = 2
        col_dim = 3
        scale=1.1
              
    # define palette colorws & subplots
    flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    fig, ax = plt.subplots(row_dim, col_dim, figsize=fsize, sharey=False, sharex=False)
    axes = ax.ravel()
    sns.set(font_scale=scale)
              
    for ax_n, y, hue, df in zip(axes, y_list, hue_list, df_list):
        # creating plot
        sns.lineplot(x='timepoint', y=y, hue=hue, data=df, legend='full', ax=ax_n,
                    palette=sns.color_palette(flatui[:len(df[hue].unique())]),)
        # manipulating axes, legend
        plt.setp(ax_n.get_xticklabels(), rotation=30)
        ax_n.legend(fontsize=12, loc='upper center')
        ax_n.set_ylim(ylim_dict[y])  
          
    # make final plot blank for chr aberr fig
    if len(y_list) > 2:       
        axes[-1].grid(False)
        axes[-1].set_facecolor('white')
        axes[-1].set_yticks([])
        axes[-1].set_xticks([])
        axes[-1].set_title('')
    plt.tight_layout()
        
    # save
    if 'telo means' in y_list and '# short telomeres' in y_list:
        plt.savefig(f'../graphs/paper figures/main figs/viz clustering groups {y_list[0]} and {y_list[1]}.png', 
                    dpi=400, bbox_inches = "tight")
    else:
        plt.savefig(f'../graphs/paper figures/main figs/viz all chr aberrations.png', 
                    dpi=400, bbox_inches = "tight")
              
              
def swap_short_telos_group_number(row):
    """
    swaps all group 1 patients to group 2 & all group 2 patients to group 1
    this is done purely for visualization reasons and has no impact whatsoever
    on interpretation of results; just swapping #s for interpretability
    """
    if row == 1:
        row = 2
    elif row == 2:
        row =1
    return row