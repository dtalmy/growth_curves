'''
Data Formatting for Host-Virus Interactions Data


Created @audralhinson
Created Aug 16th, 2019
Updated Sept 10th, 2019
'''

import numpy as np
import matplotlib as plt
from pylab import *
from scipy import *
import pandas as pd
import glob
import os
import math
import matplotlib.backends.backend_pdf
#from batch_fitting_class import *

#########################################################
################## DATA PREP ############################
#########################################################

### Find all the data ####
path = '../data/algv/' # relative paths are ok as long as you maintain the directory structure

# read in filenames according to standard naming convention
Allfiles = glob.glob(path+'*.txt') # all files
files = glob.glob(path+'[!Std]*txt') # NOT standard deviations
stddevfiles = glob.glob(path+'StdDev*') # standard deviations

# function to retreive data from txt files
def read_txt_file(file):
    filename = file
    if 'StdDev' not in str(filename): 
        HVDat = pd.read_csv(filename,names = ['time','abundance'])
    else: 
        HVDat = pd.read_csv(filename, names = ['time','abundance', 'lowlim','highlim'])
    FileSplit = re.sub('../../data/','',file)
    FileSplit = re.sub('.txt','',FileSplit)
    HVDat['description'] = FileSplit
    temp_val =  HVDat['description'].str.split('_',n=1,expand=True)
    temp_val2 =  FileSplit.split("_")
    HVDat = HVDat.dropna(axis = 0, how = 'any')
    HVDat['time'] = pd.to_numeric(HVDat['time'])
    HVDat['abundance'] = pd.to_numeric(HVDat['abundance'])
    if 'Time' in str(filename):
        HVDat['time'] = HVDat.apply(lambda x: x*24)
    if 'Log' in str(filename):
        if 'StdDev' in str(filename) or 'StdErr' in str(filename):
            HVDat['lowlim']= pow(10, (HVDat['abundance'] - HVDat['lowlim']))
            HVDat['highlim']= pow(10,(HVDat['abundance'] + HVDat['highlim']))
            HVDat['abundance'] = pow(10, HVDat['abundance'])
        else:
            HVDat['abundance'] = pow(10, HVDat['abundance'])
    HVDat['time'] = HVDat['time'].apply(lambda x: (round((x/24),1)))
    HVDat['abundance'] = HVDat['abundance'].apply(lambda x: round(x,2))
    if 'StdDev' in str(filename):
        HVDat['avglim'] = HVDat[['lowlim','highlim']].mean(axis=1)
        HVDat['stddev'] = HVDat['avglim'].apply(lambda x: round(x,2))
        HVDat = HVDat.drop(['lowlim','highlim','avglim'],axis=1)
    else:
        HVDat['stddev'] = "NA"
    HVDat['raw_data'] = temp_val[0]
    HVDat['manip'] = temp_val[1]
    HVDat['group'] = temp_val2[-1]
    HVDat = HVDat.drop('description',axis=1)
    return (HVDat)

raw_df = pd.DataFrame()

for ind,tag in enumerate(Allfiles):
    temp_df = read_txt_file(tag)
    raw_df = raw_df.append(temp_df)

raw_df = raw_df.dropna()
file_types = raw_df['raw_data'].unique().tolist()
dict_df = dict.fromkeys(file_types)

exp_set_all = []
vals_all = []
all_dat = []

############################################

print('###############################')
print('Load data')
print('###############################')
for tag in file_types:
    temp_store = raw_df[raw_df['raw_data'] == tag]
    if '../data/algv/' in tag:
        tag = tag.replace('../data/algv/','')
    print(tag,' loaded')
    dict_df[tag] = pd.pivot_table(temp_store,values='abundance', index='group',columns = 'time').T.reset_index()
    infect_df = temp_store[temp_store['group']=="Infected"]
    virus_df = temp_store[temp_store['group']=="Virus"]
    htimes = infect_df['time']
    htimes = np.asarray(infect_df[['time']]).ravel()
    habund = infect_df[['abundance']]
    habund = np.asarray(infect_df[['abundance']]).ravel()
    hstd = infect_df[['stddev']]
    hstd = np.asarray(infect_df[['stddev']]).ravel()
    vtimes = virus_df[['time']]
    vtimes = np.asarray(virus_df[['time']]).ravel()
    vabund = virus_df[['abundance']]
    vabund = np.asarray(virus_df[['abundance']]).ravel()
    vstd = virus_df[['stddev']]
    vstd = np.asarray(virus_df[['stddev']]).ravel()
    if "NA" in hstd and "NA" in vstd: 
        exp_set= {'htimes': htimes,'vtimes':vtimes,'hms':habund,'vms':vabund}
    elif "NA" in hstd: 
        exp_set= {'htimes': htimes,'vtimes':vtimes,'hms':habund,'vms':vabund,'vss':vstd}
    elif "NA" in vstd: 
        exp_set= {'htimes': htimes,'vtimes':vtimes,'hms':habund,'vms':vabund,'hss':hstd}
    else: 
        exp_set= {'htimes': htimes,'vtimes':vtimes,'hms':habund,'vms':vabund,'hss':hstd,'vss':vstd}
    # davids temporary hack just to format correctly for ODElib #
    hdat = pd.DataFrame({'abundance':habund.astype(float),'time':htimes,'sigma':habund.astype(float)/2.0})
    vdat = pd.DataFrame({'abundance':vabund.astype(float),'time':vtimes,'sigma':vabund.astype(float)/2.0})
    hdat = hdat.assign(organism='S')
    vdat = vdat.assign(organism='V')
    full_tseries = pd.concat((hdat,vdat))
    full_tseries = full_tseries.assign(paper=tag)
    all_dat.append(full_tseries)
    val = str(tag)
    exp_set_all.append(exp_set.copy())
    vals_all.append(val)
    # end davids temporary formatting hack #
