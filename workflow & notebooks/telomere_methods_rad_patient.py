import sys
import os
import pandas as pd
import numpy as np
from pandas import ExcelWriter
from pandas import ExcelFile

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import seaborn as sns

from statistics import mean 
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats

from ast import literal_eval
import more_itertools

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from xgboost import XGBRegressor
from sklearn.metrics import explained_variance_score
from sklearn.metrics import median_absolute_error


def generate_dictionary_from_TeloLength_and_Chr_aberr_Data(patharg):

#     """
#     opens raw telomere length count excel files from imageJ analyses and
#     extracts the individual mean telomere lengths to make histograms;
#     opens chromosome rearrangement frequency files and extracts data
#     both telos & chr rearrangement frequencies are stored as values to their
#     sample timepoint keys, which themselves are values to patient id# key

#     i.e the data structure is:

#     dict = {
#     patient_IDnumber = 
#     {SW#A non irrad: [telos data, chr aberr data], 
#     SW#A irrad @ 4 Gy: [telos data, chr aberr data]},

#     etc.
#     }

#     i.e:

#     all_patients_dict = {
#     '1' = {
#     'SW1A non irrad': ['telomere data', 'chr aberr data'],
#     'SW1A irrad @ 4 Gy': ['telomere data', ' chr aberr data']},

#     etc. for patients 1 - 16 (less #4 missing)
#     }

#     pass the directory where the telomere length excel files (.xlsx) are located
#     """

    all_patients_dict = {}

    for file in os.scandir(patharg):
        if file.name.endswith('.xlsx') and file.name.startswith('~$') == False:
        
            try:
                df = pd.read_excel(file)

            except:
                print('File not found..')
                return -1

            print(file.name, 'data extraction in progress..') 
#                   'it works peggy!! <3 <3 !!')
           
            if 'chr' not in file.name:
                
                telo_data = extract_and_clean_telos(df, file.name)

            else:
                continue

            file = file.name.replace('.xlsx', '').rstrip()
            data_list = []
            file_chr = ''

            
            num, num2 = capture_patient_sample_ID(file)

            if 'chr' in file:
                file_chr = file
                file = file.replace('chr','').rstrip()

            if file[num:num2] not in all_patients_dict.keys():
                all_patients_dict[file[num:num2]] = {file: []}

                if len(all_patients_dict[file[num:num2]][file]) == 0:
                    all_patients_dict[file[num:num2]][file] = data_list
                    if 'chr' not in file_chr:
                        data_list.append(telo_data)
                        data_list.sort()
                    elif 'chr' in file_chr:
                        data_list.append(chr_data)
                        data_list.sort()

                elif len(all_patients_dict[file[num:num2]][file]) == 1:
                    if 'chr' not in file_chr:
                        data_list.append(telo_data)
                        data_list.sort()
                    elif 'chr' in file_chr:
                        data_list.append(chr_data)
                        data_list.sort()

            elif file[num:num2] in all_patients_dict.keys():
                if file in all_patients_dict[file[num:num2]]:
                    if 'chr' not in file_chr:
                        all_patients_dict[file[num:num2]][file].append(telo_data)
                        all_patients_dict[file[num:num2]][file].sort()
                    elif 'chr' in file_chr:
                        all_patients_dict[file[num:num2]][file].append(chr_data)
                        all_patients_dict[file[num:num2]][file].sort()

                elif file not in all_patients_dict[file[num:num2]]:     
                    all_patients_dict[file[num:num2]][file] = data_list
                    if 'chr' not in file_chr:
                        all_patients_dict[file[num:num2]][file].append(telo_data)
                        all_patients_dict[file[num:num2]][file].sort()
                    elif 'chr' in file_chr:
                        all_patients_dict[file[num:num2]][file].append(chr_data)
                        all_patients_dict[file[num:num2]][file].sort()
                        
    print('completed file collection')
    return all_patients_dict



def generate_dataframe_from_dict_and_generate_histograms_stats(all_patients_dict, option='no graphs'):

    data = []
    print('To display graphs pass the value "yes graphs" to the function',
          'otherwise default option="no graphs"')

    for i in range(1,17):
        if str(i) in all_patients_dict.keys():
            for sample in sorted(all_patients_dict[str(i)].keys()):
                telos = all_patients_dict[str(i)][sample][0]
                # chr = all_patients_dict[str(i)[timepoint][1]
                chr_d = 'chr data'
                working_status = 'IT WORKS PEGGY <333'
            
        
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
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos, 'rsamp')
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '1 ' + 'non irrad', telos_samp.multiply(CF_mean), individ_cells, chr_d, working_status])

                elif 'irrad @ 4 Gy' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_A_irrad4Gy_name = sample
                    SW_A_irrad4Gy = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos, 'rsamp')
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '2 ' + 'irrad @ 4 Gy', telos_samp.multiply(CF_mean), individ_cells, chr_d, working_status])

                elif 'B' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_B_name = sample
                    SW_B = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos, 'rsamp')
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '3 ' + 'B', telos_samp.multiply(CF_mean), individ_cells, chr_d, working_status])
                    
                elif 'C' in sample:
                    num, num2 = capture_patient_sample_ID(sample)
                    SW_C_name = sample
                    SW_C = telos
                    telos_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, telos, 'rsamp')
                    telos_samp = telos_samp.iloc[:,0]
                    individ_cells = chunk_individual_telos_to_cells(telos_samp.multiply(CF_mean), 92)
                    data.append([sample[num:num2], '4 ' + 'C', telos_samp.multiply(CF_mean), individ_cells, chr_d, working_status])

                else:
                    print('error with making dataframe from dict..')
                    print(sample)
                    continue
            
            if option == 'yes graphs':
                
                SW_A_nonRAD_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, SW_A_nonRAD, 'rsamp')
                SW_A_irrad4Gy_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, SW_A_irrad4Gy, 'rsamp')
                SW_B_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, SW_B, 'rsamp')
                SW_C_samp = gen_missing_values_andimpute_or_randomsampledown(50, 92, SW_C, 'rsamp')

                SW_A_nonRADarray = SW_A_nonRAD_samp.to_numpy()
                SW_A_irrad4Gyarray = SW_A_irrad4Gy_samp.to_numpy()
                SW_Barray = SW_B_samp.to_numpy()
                SW_Carray = SW_C_samp.to_numpy()


                n_bins = 50
                fig, axs = plt.subplots(2, 2, sharey=True, tight_layout=False, figsize=(20, 13))

                histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, SW_A_nonRAD_samp, SW_A_nonRADarray, SW_A_nonRAD_name, 0, 0)
                histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, SW_A_irrad4Gy_samp, SW_A_nonRADarray, SW_A_irrad4Gy_name, 0, 1)
                histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, SW_B_samp, SW_A_nonRADarray, SW_B_name, 1, 0)
                histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, SW_C_samp, SW_A_nonRADarray, SW_C_name, 1, 1)

                if 'BJ1' not in sample and 'hTERT' not in sample:
                    plt.savefig(f'SW{sample[2]}_histogram.pdf')
                plt.show()
            
            else:
                continue
                
    
    all_patients_df = pd.DataFrame(data, columns=['patient id', 'timepoint', 'telo data', 'cell data', 'chr data', 'status'])
    all_patients_df['patient id'] = all_patients_df['patient id'].astype('int')
    all_patients_df = all_patients_df.sort_values(by=['patient id', 'timepoint'], ascending=True, axis=0).reset_index(drop=True)
    all_patients_df['telo means'] = all_patients_df['telo data'].apply(lambda row: np.mean(row))
    
    all_patients_df['Q1'] = 'telos nonRAD Q1 <0.25'
    all_patients_df['Q2-3'] = 'telos nonRAD Q2-3 >0.25 & <0.75'
    all_patients_df['Q4'] = 'telos nonRAD Q4 >0.75'

    return all_patients_df



def gen_missing_values_andimpute_or_randomsampledown(n_cells, telosPercell, astro_df, option=None):
    #if wanted to do for max. possible telomeres, just replace the subtraction with max telos
    # print('substracts second astro from first.. equalizing second to first')

    if astro_df.size > 4600:
        astro_dfsampled = astro_df.sample(4600)
        return astro_dfsampled

    if astro_df.size > 25 and astro_df.size <= 2300:
        missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
        rsampled = astro_df.sample(missing_data_difference, replace=True, random_state=28)
        concat_ed = pd.concat([rsampled, astro_df], sort=False)
        np.random.shuffle(concat_ed.to_numpy())
        concat_ed.reset_index(drop=True, inplace=True)
        return concat_ed

    if astro_df.size > 25 and astro_df.size < 4600:
        missing_data_difference = abs( (n_cells * telosPercell) - astro_df.size )
        if option == 'rsamp':
            rsampled = astro_df.sample(missing_data_difference, random_state=28)
            concat_ed = pd.concat([rsampled, astro_df], sort=False)
            np.random.shuffle(concat_ed.to_numpy())
            concat_ed.reset_index(drop=True, inplace=True)
            return concat_ed
        else:
            return astro_df
    else:
        return astro_df
    
    
    
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
    
    
    
def histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, astroDF, astroquartile, astroname, axsNUMone, axsNUMtwo):

        astroarray = astroDF.to_numpy()

        N, bins, patches = axs[axsNUMone,axsNUMtwo].hist(astroarray, bins=n_bins, range=(0, 400), edgecolor='black')

        for a in range(len(patches)):
#             print(bins)
#             [  0.   8.  16.  24.  32.  40.  48.  56.  64.  72.  80.  88.  96. 104.
        #  112. 120. 128. 136. 144. 152. 160. 168. 176. 184. 192. 200. 208. 216.
        #  224. 232. 240. 248. 256. 264. 272. 280. 288. 296. 304. 312. 320. 328.
        #  336. 344. 352. 360. 368. 376. 384. 392. 400.]

            if bins[a] <= np.quantile(astroquartile, 0.25):
                patches[a].set_facecolor('#fdff38')
            elif np.quantile(astroquartile, 0.25) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.50):
                patches[a].set_facecolor('#d0fefe')
            elif np.quantile(astroquartile, 0.50) < bins[a] and bins[a] <= np.quantile(astroquartile, 0.75):
                patches[a].set_facecolor('#d0fefe')
            elif bins[a] > np.quantile(astroquartile, 0.75): 
                patches[a].set_facecolor('#ffbacd')

        axs[axsNUMone,axsNUMtwo].set_title('Histogram of ' + astroname + 's Telomeres')
        axs[axsNUMone,axsNUMtwo].set_xlabel('Bins of Individ. Telomeres')
        axs[axsNUMone,axsNUMtwo].set_ylabel('Freqs of Individ. Telomeres')
        axs[axsNUMone,axsNUMtwo].xaxis.set_major_locator(plt.MaxNLocator(12))
            
        


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
    mean_values_of_individual_telomere_lengths = mean_values_of_individual_telomere_lengths.drop(labels=[5, 192, 379, 566, 753, 940, 1127, 1314,
        1501, 1688, 1875, 2062, 2249, 2436, 2623, 2810, 2997, 3184, 3371, 3558, 3745, 3932, 4119, 4306, 4493, 4680, 4867, 5054, 5241, 5428,
        5615, 5802, 5989, 6176, 6363, 6550, 6737, 6924, 7111, 7298, 7485, 7672, 7859, 8046, 8233, 8420, 8607, 8794, 8981, 9168])

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
    
    
    
### FIND QUARTILES OF NON IRRAD TIMEPOINT & MAKE BASELINE..
### find individual telomeres below the 0.25 percentile (a), between
### the 0.25 & 0.75 percentile (b), & above the 0.75 percentile (c)

def quartile_cts_rel_to_df1(df1, df2):
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    
    quartile_1 = df2[df2 <= df1.quantile(0.25)].count()
    
    quartile_2_3 = df2[(df2 > df1.quantile(0.25)) & (df2 < df1.quantile(0.75))].count()

    quartile_4 = df2[df2 >= df1.quantile(0.75)].count()
    
    return float(quartile_1.values), float(quartile_2_3.values), float(quartile_4.values)
#     return quartile_1, quartile_2_3, quartile_4



### LOOP THROUGH DATAFRAME FOR EACH PATIENT, ESTABLISH BASELINE QUARTILES FOR INDIVIDUAL TELOMERES USING NON IRRAD 
### SAMPLE TIMEPOINT.. THEN DETERMINES FOR EACH TIMEPOINT (irrad 4 Gy, B, C) HOW MANY TELOMERES REMAIN IN THOSE 
### QUARTILES... FILLS OUT Q1, Q2-3, Q4 COLUMNS..

def calculate_apply_teloQuartiles_dataframe(all_patients_df):
    
    q1_row, q2_3_row, q4_row = 7, 8, 9

    for i, row in all_patients_df.iterrows():
        if 'non irrad' in row[1]:
            nonRAD = row[2]
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, nonRAD))

        elif 'irrad @ 4 Gy' in row[1]:
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, row[2]))

        elif 'B' in row[1]:
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, row[2]))

        elif 'C' in row[1]:
            all_patients_df.iat[i, q1_row], all_patients_df.iat[i, q2_3_row], all_patients_df.iat[i, q4_row] = (quartile_cts_rel_to_df1(nonRAD, row[2]))

        else:
            print('unknown label in row[1] of the all patients df.. please check patient timepoint names')
            
    return all_patients_df



def score_model_accuracy_metrics(models, X, y):
    
    score_list = []
    
    X_train, X_valid, y_train, y_valid = train_test_split(X, y, random_state=1)
    
    for model_name in models:
    
        if model_name == 'XGBRegressor':
            model = XGBRegressor(objective='reg:squarederror', random_state=0)
        elif model_name == 'RandomForestRegressor':
            model = RandomForestRegressor(n_estimators=100, random_state=1)
            
        model.fit(X_train, y_train)
        predict_y = model.predict(X_valid)
        mae = mean_absolute_error(y_valid, predict_y)
        evs = explained_variance_score(y_valid, predict_y)
        score_list.append([model, model_name, mae, evs])
        
    score_df = pd.DataFrame(score_list, columns=['model', 'model name', 'Mean Absolute Error', 'Explained Variance'])
    return score_df, score_list


def histogram_plot_groups(x=None, data=None, 
                                 groupby=None, iterable=None):
    
    group_df = data.groupby(groupby)
    
    for item in iterable:
        plot_df = group_df.get_group(item)
        
        non_irrad = plot_df[plot_df['timepoint'] == '1 non irrad'][x]
        irrad_4_Gy = plot_df[plot_df['timepoint'] == '2 irrad @ 4 Gy'][x]
        three_B = plot_df[plot_df['timepoint'] == '3 B'][x]
        four_C = plot_df[plot_df['timepoint'] == '4 C'][x]

        n_bins = 70
        fig, axs = plt.subplots(2, 2, sharey=True, tight_layout=False, figsize=(20, 13))
        
        ax = sns.set_style(style="darkgrid",rc= {'patch.edgecolor': 'black'})
        ax = sns.set(font_scale=1)
        
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, non_irrad, non_irrad, f'patient #{item} 1 non rad', 0, 0)
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, irrad_4_Gy, non_irrad, f'patient #{item} 2 irrad @ 4 Gy', 0, 1)
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, three_B,  non_irrad, f'patient #{item} 3 B', 1, 0)
        telo_mrp.histogram_stylizer_divyBins_byQuartile(fig, axs, n_bins, four_C,  non_irrad, f'patient #{item} 4 C', 1, 1)
        
        

def color_seaborn_histogram(data, ax, bins):
    
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
        if bin_vals[a] < np.quantile(test, 0.25):
            ax.patches[a].set_facecolor('#fdff38')
        elif np.quantile(test, 0.25) < bin_vals[a] and bin_vals[a] <= np.quantile(test, 0.50):
            ax.patches[a].set_facecolor('#d0fefe')
        elif np.quantile(test, 0.50) < bin_vals[a] and bin_vals[a] <= np.quantile(test, 0.75):
            ax.patches[a].set_facecolor('#d0fefe')
        elif bin_vals[a] > np.quantile(test, 0.75): 
            ax.patches[a].set_facecolor('#ffbacd')
            
            
            
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
        'error... unknown row id'
        print(row)

        
        
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