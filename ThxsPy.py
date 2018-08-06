# -*- coding: utf-8 -*-
'''
This script reads ICP-MS output and calculates U/Th activities, which can then be plugged into the 'master' spreadsheet
Yuxin Zhou
yzhou@ldeo.columbia.edu
Next on todo list:
    find outlier in data and give warning to user
'''
import numpy as np
import numpy.ma as ma
from scipy import stats # for linear regression
from Tkinter import Tk
from tkFileDialog import askopenfilenames, asksaveasfilename
import tkMessageBox
import sys
import pandas as pd
import ctypes
import platform
#import subprocess

# check OS
if platform.system() == 'Windows':
    spike_answer = ctypes.cdll.user32.MessageBoxA(0, "Are you using 2006-2 UTh spike? If not, click no and search \'unspike\' in script and change its values", "2006-2 spike?", 4)
    if spike_answer == 7:
        sys.exit()
elif platform.system() == 'Darwin':
    window = Tk()
    window.wm_withdraw()
    tkMessageBox.showinfo(title="2006-2 spike?", message="Are you using 2006-2 UTh spike? If not, click no and search \'unspike\' in script and change its values")

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
file_names = askopenfilenames(title="Select all the ICPMS output files and a \'sample_info' file") # show an "Open" dialog box and return the path to the selected file

def return_five_point_avg(file_name):
    # start reading from row 12, which are name/unit/blank
    txt_handle = np.genfromtxt(file_name, delimiter='\t', skip_header=12)
    # get rid of first column (mass) and last column (nan)
    txt_handle = txt_handle[:,1:-1]
    # If not blank, check and remove outliers
    if 'Blank' not in file_name and 'blank' not in file_name:
        txt_handle = reject_outliers(txt_handle)
    # average accros five points
    five_point_avg = ma.mean(txt_handle.reshape(len(txt_handle)/5, 5, -1),axis=1)
    # A second check for outliers after the five point average, except when the file is Blank
    if 'Blank' not in file_name and 'blank' not in file_name:
        print file_name + " # outliers: " + str(ma.count_masked(five_point_avg))
        return reject_outliers(five_point_avg)
    else:
        return five_point_avg
    
def reject_outliers(data, m = 3.):
    # from https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
#    # median approach
#    d = np.abs(data - np.median(data, axis=1)[:,None])
#    mdev = np.median(d,axis=1)
#    s = d/mdev[:,None]# if mdev!=0 else 0.
    # avg approach
    d = np.abs(data - np.mean(data,axis=1)[:,None])
    s = d/np.mean(data,axis=1)[:,None]
    return ma.array(data,mask=np.greater(s,m))

#%% process blanks and stds. Calculate tailcrxn slope and intercept
names = [name for name in file_names if 'Blank' in name or 'blank' in name or
     'Th_std' in name or 'U_std' in name]
if not names:
    raise RuntimeError('No blank or std files found!')

# set up lists for tail corrections
# the three lists are for 238, 236, 234
U_std_tailCrxn = [[],[],[]]
blank_U_tailCrxn = [[],[],[]]
# the four lists are for 232, 229, 230, 234
Th_std_tailCrxn = [[],[],[],[]]
blank_Th_tailCrxn = [[],[],[],[]]

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)

    two_hundred_run_avg = np.mean(five_point_avg, axis=1)
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
Th234_TailCrxn = np.concatenate((Th_std_tailCrxn[3], blank_Th_tailCrxn[3]))
slopes_tailCrxn[3], intercepts_tailCrxn[3], correlations_tailCrxn[3] = stats.linregress(Th232_TailCrxn, Th234_TailCrxn)[:3]

#%% SRM_a
names = [name for name in file_names if 'SRM' in name and 'analog' in name]
if not names:
    raise RuntimeError('No SRM_a files found!')

# set up lists to store the SRM_a
SRM_a_238235_avg = []
SRM_a_238235_std = []
SRM_a_238235_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238235_avg = np.mean(five_point_avg[2
                                                    ,:]/five_point_avg[1,:])
    SRM_a_238235_avg.append(two_hundred_run_238235_avg)
    two_hundred_run_238235_std = np.std(five_point_avg[2
                                                   ,:]/five_point_avg[1,:])/np.sqrt(200)
    SRM_a_238235_std.append(two_hundred_run_238235_std)
    two_hundred_run_238235_RSD = two_hundred_run_238235_std/two_hundred_run_238235_avg
    SRM_a_238235_RSD.append(two_hundred_run_238235_RSD)
    
#%% SRM_c
names = [name for name in file_names if 'SRM' in name and 'massbias' in name]
if not names:
    raise RuntimeError('No SRM massbias files found!')
    
# set up lists to store the SRM_c
SRM_c_238235_avg = []
SRM_c_238235_std = []
SRM_c_238235_RSD = []

for file_name in names:
    five_point_avg = return_five_point_avg(file_name)
    
    two_hundred_run_238235_avg = np.mean(five_point_avg[2
                                                    ,:]/five_point_avg[1,:])
    SRM_c_238235_avg.append(two_hundred_run_238235_avg)
    two_hundred_run_238235_std = np.std(five_point_avg[2
                                                   ,:]/five_point_avg[1,:])/np.sqrt(200)
    SRM_c_238235_std.append(two_hundred_run_238235_std)
    two_hundred_run_238235_RSD = two_hundred_run_238235_std/two_hundred_run_238235_avg
    SRM_c_238235_RSD.append(two_hundred_run_238235_RSD)

#%% sample results
# set up the 2d array as in export spreadsheet
# Columns: 238/236_avg	238/236_RSD	235/236_avg	235/236_RSD	234/236_avg	234/236_RSD	230/229_avg	230/229_stdev	232/229_avg	232/229_stdev
# Rows: UTh1-20
export = np.zeros((20,10))

# if this is UTh data file
names = [name for name in file_names if 'UTh.txt' in name]
if not names:
    raise RuntimeError('No UTh files found!')
names.sort()
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
    two_hundred_run_238236_avg = np.mean(five_point_avg[-1
                                                        ,:]/five_point_avg[-2,:])
    export[i,0] = two_hundred_run_238236_avg
    two_hundred_run_238236_std = np.std(five_point_avg[-1
                                                        ,:]/five_point_avg[-2,:])/np.sqrt(200)
    two_hundred_run_238236_RSD = two_hundred_run_238236_std/two_hundred_run_238236_avg
    export[i,1] = two_hundred_run_238236_RSD
    # 235/236 U
    two_hundred_run_235236_avg = np.mean(five_point_avg[-3
                                                        ,:]/five_point_avg[-2,:])
    export[i,2] = two_hundred_run_235236_avg
    two_hundred_run_235236_std = np.std(five_point_avg[-3
                                                        ,:]/five_point_avg[-2,:])/np.sqrt(200)
    two_hundred_run_235236_RSD = two_hundred_run_235236_std/two_hundred_run_235236_avg
    export[i,3] = two_hundred_run_235236_RSD
    # 234/236 U
    two_hundred_run_234236_avg = np.mean(five_point_avg[3
                                                        ,:]/five_point_avg[-2,:])
    export[i,4] = two_hundred_run_234236_avg
    two_hundred_run_234236_std = np.std(five_point_avg[3
                                                        ,:]/five_point_avg[-2,:])/np.sqrt(200)
    two_hundred_run_234236_RSD = two_hundred_run_234236_std/two_hundred_run_234236_avg
    export[i,5] = two_hundred_run_234236_RSD
    # 230/229 Th
    two_hundred_run_230229_avg = np.mean(five_point_avg[1
                                                        ,:]/five_point_avg[0,:])
    export[i,6] = two_hundred_run_230229_avg
    two_hundred_run_230229_std = np.std(five_point_avg[1
                                                        ,:]/five_point_avg[0,:])/np.sqrt(200)
    two_hundred_run_230229_RSD = two_hundred_run_230229_std/two_hundred_run_230229_avg
    export[i,7] = two_hundred_run_230229_RSD
    # 232/229 Th
    two_hundred_run_232229_avg = np.mean(five_point_avg[2
                                                        ,:]/five_point_avg[0,:])
    export[i,8] = two_hundred_run_232229_avg
    two_hundred_run_232229_std = np.std(five_point_avg[2
                                                        ,:]/five_point_avg[0,:])/np.sqrt(200)
    two_hundred_run_232229_RSD = two_hundred_run_232229_std/two_hundred_run_232229_avg
    export[i,9] = two_hundred_run_232229_RSD
    
#%% todo ez reduction

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
SRM_c = np.mean(SRM_c_238235_avg)
SRM_c_RSD = np.sqrt((np.sum((np.array(SRM_c_238235_avg) * np.array(SRM_c_238235_RSD))**2)))/3/SRM_c
accepted_238235 = 137.55
accepted_238235_RSD = 0.50*0.01
mass_bias_per_amu = (SRM_c/accepted_238235-1)/3
mass_bias_per_amu_RSD = np.sqrt((SRM_c_RSD**2+accepted_238235_RSD**2))

# analog counting gain
SRM_a = np.mean(SRM_a_238235_avg)
SRM_a_RSD = np.sqrt((np.sum((np.array(SRM_a_238235_avg) * np.array(SRM_a_238235_RSD))**2)))/3/SRM_a
counting_gain = SRM_a/SRM_c
counting_gain_RSD = np.sqrt(SRM_a_RSD**2+SRM_c_RSD**2)

##Magic
# mass bias correction
mass_difference = np.array([2,-1,-2,1,3])
mass_bias_crxn = 1+mass_difference*mass_bias_per_amu
mass_bias_crxn_RSD = mass_bias_per_amu_RSD*abs(mass_difference)
for i in range(5):
    export[:,i*2] *= mass_bias_crxn[i]
    export[:,i*2+1] = np.sqrt(export[:,i*2+1]**2+mass_bias_crxn_RSD[i]**2)

# counting gain crxn
for i in [0,8]:
    export[:,i] /= counting_gain
    export[:,i+1] = np.sqrt(export[:,i+1]**2+counting_gain_RSD**2)
#%%
# unspike
unspike_matrix = np.array([[1.535934579,0.004656157,0.0030],[2.511807533,0.005569552,0.0022]])
for i,weight in [(0,238),(1,235),(2,234)]:
	if sample_info_type == 'txt':
		export[:,i*2] = sample_info['f3']/1000*unspike_matrix[1,0]*export[:,i*2]*weight/236
	elif sample_info_type == 'xlsx':
		export[:,i*2] = sample_info[3]/1000*unspike_matrix[1,0]*export[:,i*2]*weight/236
	export[:,i*2+1] = np.sqrt(unspike_matrix[1,2]**2+export[:,i*2+1]**2)
for i,weight in [(3,230),(4,232)]:
	if sample_info_type == 'txt':
		export[:,i*2] = sample_info['f3']/1000*unspike_matrix[0,0]*export[:,i*2]*weight/229
	elif sample_info_type == 'xlsx':
		export[:,i*2] = sample_info[3]/1000*unspike_matrix[0,0]*export[:,i*2]*weight/229
	export[:,i*2+1] = np.sqrt(unspike_matrix[0,2]**2+export[:,i*2+1]**2)
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
	blank_index = np.argwhere(sample_info[0]=='BLANK')
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
export_data_df = pd.DataFrame(data=export,index=np.arange(20),columns=['238U dpm/g',	'238U dpm/g 2 sigma',	'234U dpm/g',	'234U dpm/g 2 sigma',	'230Th dpm/g',	'230Th dpm/g 2 sigma',	'232Th dpm/g',	'232Th dpm/g 2 sigma'])
if sample_info_type == 'txt':
	sample_name_df = pd.DataFrame({'Sample name':sample_info['f0']})
elif sample_info_type == 'xlsx':
	sample_name_df = pd.DataFrame({'Sample name':sample_info[0]})
export_df = pd.concat([sample_name_df,export_data_df],axis=1)

#%% save to csv
output_file_name = asksaveasfilename(title='Save the output file as')
if 'xlsx' not in output_file_name:
    output_file_name = output_file_name + '.xlsx'
export_df.to_excel(output_file_name)