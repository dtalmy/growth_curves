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
from batch_fitting_class import *

print('testing, are we collaborators?')
#########################################################
################## DATA PREP ############################
#########################################################

### Find all the data ####
path = '../data/algv/' # relative paths are ok as long as you maintain the directory structure

# read in filenames according to standard naming convention
Allfiles = glob.glob(path+'*.txt') # all files
files = glob.glob(path+'[!Std]*txt') # NOT standard deviations
stddevfiles = glob.glob(path+'StdDev*') # standard deviations

# read in the data
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
    val = str(tag)
    exp_set_all.append(exp_set.copy())
    vals_all.append(val)

#########################################################
################## THE MODELS ###########################
#########################################################

###Create the result graph pdf ####
pdf = matplotlib.backends.backend_pdf.PdfPages("../figures/ModelTester.pdf")

# define which processes to include in the model
nstates = arange(0,3) # number of infection states
chronic = False # budding of cells from infected state
suicide = False # programmed cell death by infected hosts
muinfec = False # growth by infected hosts
recover = False # allow infected hosts to recover from infection
pgrowth = False # allow the predator to grow in the absence of prey
hnutlim = False # allow nutrient limitation of the host
pnutlim = False # allow nutrient limitation of the predator

# put all the defining features into a single list comprehension
mods = [[n,chronic,suicide,muinfec,recover,pgrowth,hnutlim,pnutlim] for n in nstates]

# define labels to aid in interpretation of output
labs = ['#I states ='+str(n) for n in nstates]
labloc = np.arange(len(labs))  

#### Run the models ####

print(' ')
print('###############################')
print('Do fitting')
print('###############################')

for (inf,tag) in zip(exp_set_all[:1],vals_all[:1]):
    phi_all,beta_all=r_[[]],r_[[]]
    allmods = []
    for (mod,lab) in zip(mods,labs):
        pmod = all_mods(inf,mod, nits=1000,pits=100,burnin=500)
        print(' ')
        print('dataset: ',tag)
        print('model label: ', lab)
        print('model params: ', pmod.pnames)
        print(' ')
        if mod == mods[0]:
            figs,axes = pmod.gen_figs(tag)
        pmod.modellabel = lab
        pmod.do_fitting(axes[0])
        print('AIC: ', pmod.AIC)
        print('Adjusted R squared: ', pmod.adj_rsquared)
        print('final error: ', pmod.chi)
        print(' ')
        print('###############################')
        axes[1].plot(pmod.iterations,pmod.likelihoods)
        allmods.append(pmod)
        phi_all = append(phi_all,exp(pmod.pall[1]))
        for a in axes[0]:
            a.semilogy()
    leg = axes[0][0].legend()
    leg.draw_frame(False)
    axes[2][0].scatter(arange(len(mods)),r_[[mod.adj_rsquared for mod in allmods]])
    axes[2][1].scatter(arange(len(mods)),r_[[mod.AIC for mod in allmods]])
    for a in axes[2]:
        a.set_xticks(labloc)
        a.set_xticklabels(labs, fontweight = "bold")
        a.tick_params(direction = "in")
    for f in figs:
        pdf.savefig(f)
        close(f)
pdf.close()

