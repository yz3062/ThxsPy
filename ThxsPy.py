# -*- coding: utf-8 -*-
'''
This script reads ICP-MS output and calculates U/Th activities, which can then be plugged into the 'master' spreadsheet
Yuxin Zhou
yzhou@ldeo.columbia.edu
'''
# If directly import mpl, crashes. Source: https://stackoverflow.com/questions/32019556/matplotlib-crashing-tkinter-application/34109240#34109240
import matplotlib
# matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy import stats # for linear regression
import tkinter as tk
from tkinter import filedialog
import sys
import pandas as pd
import pprint

#import subprocess

spike_answer = str(input("Are you using 2006-2 UTh spike? If not, click no and search \'unspike\' in script and change its values. [y] or n:") or 'y')
if spike_answer == 'n':
    sys.exit()
figure_answer = str(input("Do you want to inspect ICPMS raw output in figures?[y] or n:") or 'y')

root = tk.Tk()
root.withdraw() # we don't want a full GUI, so keep the root window from appearing
file_names = filedialog.askopenfilenames(title="Select all the ICPMS output files and a \'sample_info' file") # show an "Open" dialog box and return the path to the selected file

def return_five_point_avg(file_name):
    # read txt as csv, using tab as separator
    txt_handle = pd.read_csv(file_name,sep='\t',header=None)
    txt_handle.dropna(how='all',axis='index',inplace=True) # drop the rows where all elements are missing
    txt_handle.dropna(how='all',axis='columns',inplace=True) # drop the columns where all elements are missing
    txt_handle.reset_index(drop=True,inplace=True) # index at this point doesn't start with 0,1,2 etc because some rows were deleted. This will reset index, and drop=True will prevent a new column named "index" be created
    txt_handle.drop([0,1,2],inplace=True) # remove the first three rows 
    txt_handle.reset_index(drop=True,inplace=True)
    txt_handle = txt_handle.astype(float)
    if figure_answer == 'y':
        txt_handle_r = txt_handle.transpose() # create a transposed version
        txt_handle_r.columns = txt_handle_r.iloc[0] # name the columns with the first row (mass)
        txt_handle_r.drop(txt_handle_r.index[0],inplace=True) # drop the first row (mass)
        txt_handle_r.plot(sharex=True,title=file_name)
        plt.savefig(file_name+'.png')
    txt_handle.set_index(txt_handle[0],inplace=True) # set index as mass
    txt_handle.drop(columns=[0],inplace=True) # drop the mass column
    txt_handle = reject_outliers(txt_handle)
    # average accros multiple masses of the same element
    txt_handle.set_index(np.floor(txt_handle.index),inplace=True)
    five_point_avg = txt_handle.groupby(txt_handle.index).mean()
    five_point_avg_2nd_outlier_detection = reject_outliers(five_point_avg)
    masekd_array = ma.masked_invalid(five_point_avg_2nd_outlier_detection.values)
    print(file_name + ' # outliers: ' + str(np.count_nonzero(~np.isnan(masekd_array))))
    return masekd_array
    
def reject_outliers(data, m = 2.):
    '''
    from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
    data is expected to be a pandas dataframe, where each row is a number of measurements on the same mass
    '''
    d = data.subtract(data.median(axis=1),axis='index').abs()
    mdev = d.median(axis=1)
    s = d.divide(mdev,axis='index')
    return data.mask(s>m)

#%% process blanks and stds. Calculate tailcrxn slope and intercept
names = [name for name in file_names if ('Blank' in name or 'blank' in name or
     'Th_std' in name or 'U_std' in name) and 'SRM' not in name]
if not names:
    raise RuntimeError('No blank or std files found!')
print("Identified the following files as either blank or U_std or Th_std:")
pprint.pprint(names)
print('\n')
# set up lists for tail corrections
# the three lists are for 238, 236, 234
U_std_tailCrxn = [[],[],[]]
blank_U_tailCrxn = [[],[],[]]
# the four lists are for 232, 229, 230, 234
Th_std_tailCrxn = [[],[],[],[]]
blank_Th_tailCrxn = [[],[],[],[]]

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)

    two_hundred_run_avg = ma.mean(five_point_avg, axis=1)
    if 'Blank' in file_name or 'blank' in file_name:
        blank_U_tailCrxn[0].append(two_hundred_run_avg[-1]) # U238
        blank_U_tailCrxn[1].append(two_hundred_run_avg[-2]) # U236
        blank_U_tailCrxn[2].append(two_hundred_run_avg[3]) # U234
        blank_Th_tailCrxn[0].append(two_hundred_run_avg[2]) # Th232
        blank_Th_tailCrxn[1].append(two_hundred_run_avg[0]) # Th229
        blank_Th_tailCrxn[2].append(two_hundred_run_avg[1]) # Th230
        blank_Th_tailCrxn[3].append(two_hundred_run_avg[3]) # Th234
    if 'U_std' in file_name:
        U_std_tailCrxn[0].append(two_hundred_run_avg[-1]) # U238
        U_std_tailCrxn[1].append(two_hundred_run_avg[-2]) # U236
        U_std_tailCrxn[2].append(two_hundred_run_avg[3]) # U234
    if 'Th_std' in file_name:
        Th_std_tailCrxn[0].append(two_hundred_run_avg[2]) # Th232
        Th_std_tailCrxn[1].append(two_hundred_run_avg[0]) # Th229
        Th_std_tailCrxn[2].append(two_hundred_run_avg[1]) # Th230
        Th_std_tailCrxn[3].append(two_hundred_run_avg[3]) # Th234

# now do TailCrxn before processing UTh data file
# two arrays to store intercepts and slopes in the sequence of
# 236 vs 238, 234 vs 238, 229 vs 232, 230 vs 232, 234 vs 232
intercepts_tailCrxn = np.zeros(4)
slopes_tailCrxn = np.zeros(4)
correlations_tailCrxn = np.zeros(4)
U238_tailCrxn = np.concatenate((U_std_tailCrxn[0], blank_U_tailCrxn[0]))
U236_tailCrxn = np.concatenate((U_std_tailCrxn[1], blank_U_tailCrxn[1]))
slopes_tailCrxn[0], intercepts_tailCrxn[0], correlations_tailCrxn[0] = stats.linregress(U238_tailCrxn, U236_tailCrxn)[:3]
#U234_tailCrxn = np.concatenate((U_std_tailCrxn[2], blank_U_tailCrxn[2]))
#slopes_tailCrxn[1], intercepts_tailCrxn[1], correlations_tailCrxn[1] = stats.linregress(U238_tailCrxn, U234_tailCrxn)[:3]
Th232_TailCrxn = np.concatenate((Th_std_tailCrxn[0], blank_Th_tailCrxn[0]))
Th229_TailCrxn = np.concatenate((Th_std_tailCrxn[1], blank_Th_tailCrxn[1]))
slopes_tailCrxn[1], intercepts_tailCrxn[1], correlations_tailCrxn[1] = stats.linregress(Th232_TailCrxn, Th229_TailCrxn)[:3]
Th230_TailCrxn = np.concatenate((Th_std_tailCrxn[2], blank_Th_tailCrxn[2]))
slopes_tailCrxn[2], intercepts_tailCrxn[2], correlations_tailCrxn[2] = stats.linregress(Th232_TailCrxn, Th230_TailCrxn)[:3]
U234_TailCrxn = np.concatenate((Th_std_tailCrxn[3], blank_Th_tailCrxn[3]))
slopes_tailCrxn[3], intercepts_tailCrxn[3], correlations_tailCrxn[3] = stats.linregress(Th232_TailCrxn, U234_TailCrxn)[:3]

#%% SRM_a blank
names = [name for name in file_names if 'SRM' in name and 'analog' in name and ('blank' in name or 'Blank' in name or 'BLANK' in name)]
if not names:
    SRM_a_blank_flag = False
else:
    SRM_a_blank_flag = True
    print("Identified the following files as SRM_a blanks:")
    pprint.pprint(names)
    print('\n')
    
# set up lists to store the 3 SRM_a_blanks
SRM_a_238_blank_avg = []
SRM_a_235_blank_avg = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238_avg = ma.mean(five_point_avg[2,:])
    two_hundred_run_235_avg = ma.mean(five_point_avg[1,:])
    SRM_a_238_blank_avg.append(two_hundred_run_238_avg)
    SRM_a_235_blank_avg.append(two_hundred_run_235_avg)

#%% SRM_a
names = [name for name in file_names if 'SRM' in name and 'analog' in name and 'blank' not in name and 'Blank' not in name and 'BLANK' not in name]
if not names:
    raise RuntimeError('No SRM_a files found!')
print("Identified the following files as SRM_a:")
pprint.pprint(names)
print('\n')

# set up lists to store the SRM_a
SRM_a_238_avg = []
SRM_a_235_avg = []
SRM_a_238235_avg = []
SRM_a_238235_std = []
SRM_a_238235_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_SRM_a_238_avg = ma.mean(five_point_avg[2,:])
    two_hundred_run_SRM_a_235_avg = ma.mean(five_point_avg[1,:])
    SRM_a_238_avg.append(two_hundred_run_SRM_a_238_avg)
    SRM_a_235_avg.append(two_hundred_run_SRM_a_235_avg)
    two_hundred_run_238235_avg = ma.mean(five_point_avg[2
                                                    ,:]/five_point_avg[1,:])
    two_hundred_run_238235_std = ma.std(five_point_avg[2
                                                   ,:]/five_point_avg[1,:])/ma.sqrt(five_point_avg.shape[1])
    SRM_a_238235_std.append(two_hundred_run_238235_std)
    two_hundred_run_238235_RSD = two_hundred_run_238235_std/two_hundred_run_238235_avg
    SRM_a_238235_RSD.append(two_hundred_run_238235_RSD)
if SRM_a_blank_flag:
    SRM_a_238235_avg = (SRM_a_238_avg - ma.mean(SRM_a_238_blank_avg)) / (SRM_a_235_avg - ma.mean(SRM_a_235_blank_avg))
else:
    SRM_a_238235_avg = ma.array(SRM_a_238_avg) / ma.array(SRM_a_235_avg)
#%% SRM_c blank
names = [name for name in file_names if 'SRM' in name and 'massbias' in name and ('blank' in name or 'Blank' in name or 'BLANK' in name)]
if not names:
    SRM_c_blank_flag = False
else:
    SRM_c_blank_flag = True
    print("Identified the following files as SRM_c blanks:")
    pprint.pprint(names)
    print('\n')
    
# set up lists to store the 3 SRM_a_blanks
SRM_c_238_blank_avg = []
SRM_c_235_blank_avg = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238_avg = ma.mean(five_point_avg[2,:])
    two_hundred_run_235_avg = ma.mean(five_point_avg[1,:])
    SRM_c_238_blank_avg.append(two_hundred_run_238_avg)
    SRM_c_235_blank_avg.append(two_hundred_run_235_avg)
    
#%% SRM_c
names = [name for name in file_names if 'SRM' in name and 'massbias' in name and 'blank' not in name and 'Blank' not in name and 'BLANK' not in name]
if not names:
    raise RuntimeError('No SRM massbias files found!')
print("Identified the following files as either SRM_c:")
pprint.pprint(names)
print('\n')    

# set up lists to store the SRM_c
SRM_c_238_avg = []
SRM_c_235_avg = []
SRM_c_238235_avg = []
SRM_c_238235_std = []
SRM_c_238235_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_SRM_c_238_avg = ma.mean(five_point_avg[2,:])
    two_hundred_run_SRM_c_235_avg = ma.mean(five_point_avg[1,:])
    SRM_c_238_avg.append(two_hundred_run_SRM_c_238_avg)
    SRM_c_235_avg.append(two_hundred_run_SRM_c_235_avg)
    two_hundred_run_238235_avg = ma.mean(five_point_avg[2
                                                    ,:]/five_point_avg[1,:])
    two_hundred_run_238235_std = ma.std(five_point_avg[2
                                                   ,:]/five_point_avg[1,:])/ma.sqrt(five_point_avg.shape[1])
    SRM_c_238235_std.append(two_hundred_run_238235_std)
    two_hundred_run_238235_RSD = two_hundred_run_238235_std/two_hundred_run_238235_avg
    SRM_c_238235_RSD.append(two_hundred_run_238235_RSD)
if SRM_c_blank_flag:
    SRM_c_238235_avg = (SRM_c_238_avg - ma.mean(SRM_c_238_blank_avg)) / (SRM_c_235_avg - ma.mean(SRM_c_235_blank_avg))
else:
    SRM_c_238235_avg = ma.array(SRM_c_238_avg) / ma.array(SRM_c_235_avg)
#%% sample results

# if this is UTh data file
names = [name for name in file_names if 'UTh.txt' in name]
if not names:
    raise RuntimeError('No UTh files found!')
names.sort()
print("Identified the following files as sample files:")
pprint.pprint(names)
print('\n')
# set up the 2d array as in export spreadsheet
# Columns: 238/236_avg    238/236_RSD    235/236_avg    235/236_RSD    234/236_avg    234/236_RSD    230/229_avg    230/229_stdev    232/229_avg    232/229_stdev
# Rows: UTh1-num_samples
num_sample = len(names)
export = np.zeros((num_sample,10))

for i, file_name in enumerate(names):
    five_point_avg = return_five_point_avg(file_name)
    
    # first correct for tailing
    # correct 229
    five_point_avg[0,:] -= slopes_tailCrxn[1] * five_point_avg[2,:] + intercepts_tailCrxn[1]
    # correct 230
    five_point_avg[1,:] -= slopes_tailCrxn[2] * five_point_avg[2,:] + intercepts_tailCrxn[2]
    # correct 234, if negative set to 0
    five_point_avg[3,:] -= slopes_tailCrxn[3] * five_point_avg[2,:] + intercepts_tailCrxn[3]
    five_point_avg[3,:][five_point_avg[3,:] < 0] = 0
    # correct 236
    five_point_avg[-2,:] -= slopes_tailCrxn[0] * five_point_avg[-1,:] + intercepts_tailCrxn[0]

    # calcualte the ratios
    # 238/236 U
    two_hundred_run_238236_avg = ma.mean(five_point_avg[-1
                                                        ,:]/five_point_avg[-2,:])
    export[i,0] = two_hundred_run_238236_avg
    two_hundred_run_238236_std = ma.std(five_point_avg[-1
                                                        ,:]/five_point_avg[-2,:])/ma.sqrt(five_point_avg.shape[1])
    two_hundred_run_238236_RSD = two_hundred_run_238236_std/two_hundred_run_238236_avg
    export[i,1] = two_hundred_run_238236_RSD
    # 235/236 U
    two_hundred_run_235236_avg = ma.mean(five_point_avg[-3
                                                        ,:]/five_point_avg[-2,:])
    export[i,2] = two_hundred_run_235236_avg
    two_hundred_run_235236_std = ma.std(five_point_avg[-3
                                                        ,:]/five_point_avg[-2,:])/ma.sqrt(five_point_avg.shape[1])
    two_hundred_run_235236_RSD = two_hundred_run_235236_std/two_hundred_run_235236_avg
    export[i,3] = two_hundred_run_235236_RSD

    # 234/236 U
    two_hundred_run_234236_avg = ma.mean(five_point_avg[3
                                                        ,:]/five_point_avg[-2,:])
    export[i,4] = two_hundred_run_234236_avg
    two_hundred_run_234236_std = ma.std(five_point_avg[3
                                                        ,:]/five_point_avg[-2,:])/ma.sqrt(five_point_avg.shape[1])
    two_hundred_run_234236_RSD = two_hundred_run_234236_std/two_hundred_run_234236_avg
    export[i,5] = two_hundred_run_234236_RSD
    # 230/229 Th
    two_hundred_run_230229_avg = ma.mean(five_point_avg[1
                                                        ,:]/five_point_avg[0,:])
    export[i,6] = two_hundred_run_230229_avg
    two_hundred_run_230229_std = ma.std(five_point_avg[1
                                                        ,:]/five_point_avg[0,:])/ma.sqrt(five_point_avg.shape[1])
    two_hundred_run_230229_RSD = two_hundred_run_230229_std/two_hundred_run_230229_avg
    export[i,7] = two_hundred_run_230229_RSD
    # 232/229 Th
    two_hundred_run_232229_avg = ma.mean(five_point_avg[2
                                                        ,:]/five_point_avg[0,:])
    export[i,8] = two_hundred_run_232229_avg
    two_hundred_run_232229_std = ma.std(five_point_avg[2
                                                        ,:]/five_point_avg[0,:])/ma.sqrt(five_point_avg.shape[1])
    two_hundred_run_232229_RSD = two_hundred_run_232229_std/two_hundred_run_232229_avg
    export[i,9] = two_hundred_run_232229_RSD
    
#%% ez reduction

# sample info. Exclude $ in file name in case that file is open
names = [name for name in file_names if 'info' in name and '$' not in name]
if len(names)>1:
    raise RuntimeError('More than one sample info file')
if not names:
    raise RuntimeError('Sample info file cannot be found. The file must have \'info\' in file name')
sample_info_type = ''
if names[0][-3:] == 'txt':
    sample_info_type = 'txt'
    try:
        sample_info = np.genfromtxt(names[0], delimiter='\t',dtype=None,skip_header=1)
    except ValueError:
        raise ValueError('In reading file ' + names[0] + ', value error!')
elif names[0][-4:] == 'xlsx':
    sample_info_type = 'xlsx'
    try:
        sample_info = pd.read_excel(names[0],header=None,skiprows=1)
    except ValueError:
        raise ValueError('In reading file ' + names[0] + ', value error!')
else:
    raise ValueError(names[0] + ' is not either a txt or excel file and cannot be processed')

## MassBiasCountGain
# mass bias
SRM_c = ma.mean(SRM_c_238235_avg)
SRM_c_RSD = ma.sqrt((ma.sum((ma.array(SRM_c_238235_avg) * ma.array(SRM_c_238235_RSD))**2)))/3/SRM_c
accepted_238235 = 137.818 # Hiess et al., 2012
accepted_238235_RSD = 0.50*0.01
mass_bias_per_amu = (accepted_238235/SRM_c-1)/3
mass_bias_per_amu_RSD = ma.sqrt((SRM_c_RSD**2+accepted_238235_RSD**2))

# analog counting gain
SRM_a = ma.mean(SRM_a_238235_avg)
SRM_a_RSD = ma.sqrt((ma.sum((ma.array(SRM_a_238235_avg) * ma.array(SRM_a_238235_RSD))**2)))/3/SRM_a
counting_gain = SRM_a/SRM_c
counting_gain_RSD = ma.sqrt(SRM_a_RSD**2+SRM_c_RSD**2)

##Magic
# mass bias correction
mass_difference = np.array([2,-1,-2,1,3])
mass_bias_crxn = 1+mass_difference*mass_bias_per_amu
mass_bias_crxn_RSD = mass_bias_per_amu_RSD*abs(mass_difference)
for i in range(5):
    export[:,i*2] *= mass_bias_crxn[i]
    export[:,i*2+1] = ma.sqrt(export[:,i*2+1]**2+mass_bias_crxn_RSD[i]**2)

# counting gain crxn
for i in [0,8]:
    export[:,i] /= counting_gain
    export[:,i+1] = ma.sqrt(export[:,i+1]**2+counting_gain_RSD**2)
#%%
# unspike
unspike_matrix = np.array([[1.535934579,0.004656157,0.0030],[2.511807533,0.005569552,0.0022]])
for i,weight in [(0,238),(1,235),(2,234)]:
    if sample_info_type == 'txt':
        export[:,i*2] = sample_info['f3']/1000*unspike_matrix[1,0]*export[:,i*2]*weight/236
    elif sample_info_type == 'xlsx':
        export[:,i*2] = sample_info[3]/1000*unspike_matrix[1,0]*export[:,i*2]*weight/236
    export[:,i*2+1] = ma.sqrt(unspike_matrix[1,2]**2+export[:,i*2+1]**2)
for i,weight in [(3,230),(4,232)]:
    if sample_info_type == 'txt':
        export[:,i*2] = sample_info['f3']/1000*unspike_matrix[0,0]*export[:,i*2]*weight/229
    elif sample_info_type == 'xlsx':
        export[:,i*2] = sample_info[3]/1000*unspike_matrix[0,0]*export[:,i*2]*weight/229
    export[:,i*2+1] = ma.sqrt(unspike_matrix[0,2]**2+export[:,i*2+1]**2)
#%%
# Sed Cncn ng/g
if sample_info_type == 'txt':
    if not (sample_info['f0']=='BLANK').any():
        raise RuntimeError('Cannot determine from sample name in sample info which sample is blank. Name it BLANK')
    blank_index = np.argwhere(sample_info['f0']=='BLANK')
    multiplication_factor=[0.001,1,1000,1000,0.001]
    for i in range(5):
        export[:,i*2] = (export[:,i*2]-export[blank_index,i*2])*multiplication_factor[i]/(sample_info['f2']/1000)
elif sample_info_type == 'xlsx':
    if not (sample_info[0]=='BLANK').any():
        raise RuntimeError('Cannot determine from sample name in sample info which sample is blank. Name it BLANK')
    # blank_index = np.argwhere(sample_info[0]=='BLANK')
    # The following line is a temp fix for a pandas bug
    # https://github.com/pandas-dev/pandas/issues/35331
    blank_index = sample_info[0][sample_info[0]=='BLANK'].index[0]
    multiplication_factor=[0.001,1,1000,1000,0.001]
    for i in range(5):
        export[:,i*2] = np.squeeze(export[:,i*2]-export[blank_index,i*2])*multiplication_factor[i]/(sample_info[2]/1000)

#%%
# Sed Cncn dpm/g
sed_cncn_dpm_matrix=[0.752049334,0.013782268,0.045747747,0.242530074]
for i,column in enumerate([0,4,6,8]):
    export[:,column]=sed_cncn_dpm_matrix[i]*export[:,column]

# two sigma
for i in range(5):
    export[:,i*2+1] = export[:,i*2+1]*export[:,i*2]*2
    
# delete 235 column
export=np.delete(export,[2,3],1)

# since numpy array can't have both string and float, converting to pandas dataframe and add sample name as the first column in export
export_data_df = pd.DataFrame(data=export,index=np.arange(num_sample),columns=['238U dpm/g',    '238U dpm/g 2 sigma',    '234U dpm/g',    '234U dpm/g 2 sigma',    '230Th dpm/g',    '230Th dpm/g 2 sigma',    '232Th dpm/g',    '232Th dpm/g 2 sigma'])
if sample_info_type == 'txt':
    sample_name_df = pd.DataFrame({'Sample name':sample_info['f0']})
elif sample_info_type == 'xlsx':
    sample_name_df = pd.DataFrame({'Sample name':sample_info[0]})
export_df = pd.concat([sample_name_df,export_data_df],axis=1)

#%% save to csv
output_file_name = filedialog.asksaveasfilename(title='Save the output file as')
if 'xlsx' not in output_file_name:
    output_file_name = output_file_name + '.xlsx'
export_df.to_excel(output_file_name)
