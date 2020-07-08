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

#########################################################
################## DATA PREP ############################
#########################################################

### Find all the data ####
path = '../data/algv/' # relative paths are ok as long as you maintain the directory structure

# read in filenames according to standard naming convention
# files = glob.glob(path+'*.txt') # all files
# files = glob.glob(path+'*[!_StdDev_Control.txt]*') # NOT standard deviations 
files = set(glob.glob(path+'*Control*txt')) - set(glob.glob(path+'*StdDev*txt')) #true not standard deviations
# stddevfiles = glob.glob(path+'StdDev*') # standard deviations

# read in the data
def read_txt_file(file):
    filename = file
    print(filename)
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

for ind,tag in enumerate(files):
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
    print(tag)
    dict_df[tag] = pd.pivot_table(temp_store,values='abundance', index='group',columns = 'time').T.reset_index()
    control_df = temp_store[temp_store['group']=="Control"]
    virus_df = temp_store[temp_store['group']=="Virus"]
    htimes = control_df['time']
    htimes = np.asarray(control_df[['time']]).ravel()
    habund = control_df[['abundance']]
    habund = np.asarray(control_df[['abundance']]).ravel()
    hstd = control_df[['stddev']]
    hstd = np.asarray(control_df[['stddev']]).ravel()
    vtimes = virus_df[['time']]
    vtimes = np.asarray(virus_df[['time']]).ravel()
    vabund = virus_df[['abundance']]
    vabund = np.asarray(virus_df[['abundance']]).ravel()
    vstd = virus_df[['stddev']]
    vstd = np.asarray(virus_df[['stddev']]).ravel()
    
    if "NA" in hstd and "NA" in vstd: 
        exp_set= {'htimes': htimes,'hms':habund} #,'vms':vabund ,'vtimes':vtimes
    elif "NA" in hstd: 
        exp_set= {'htimes': htimes,'hms':habund} #,'vms':vabund,'vss':vstd,'vtimes':vtimes
    elif "NA" in vstd: 
        exp_set= {'htimes': htimes,'hms':habund} #,'vms':vabund,'hss':hstd,'vtimes':vtimes
    else: 
        exp_set= {'htimes': htimes,'hms':habund} #,'vms':vabund,'hss':hstd,'vss':vstd,'vtimes':vtimes
    val = str(tag)
    exp_set_all.append(exp_set.copy())
    vals_all.append(val)


#########################################################
################## THE MODELS ###########################
#########################################################

###Create the result graph pdf ####
pdf = matplotlib.backends.backend_pdf.PdfPages("../figures/ModelTester.pdf")

# define the models to be tested
M1 = ['mum','phi', 'beta', 'lambd'] #Lotka Volterra
M2 = M1 + ['deltv'] 
M3 = M2 + ['alp'] #One infection class, including chronic release/budding
M4 = M2 + ['psi'] #One infection class, including chronic release/budding and host suicide
M5 = M2 + ['tau'] #Two infection classes with the same lysis and burst size
M6 = M5 + ['betal'] #Two infection classes with independent burst sizes but the same lysis rate
M7 = M6 + ['lambdl'] #Two infection classes with independent burst sizes and lysis rates
M8 = M5 + ['alp'] #Two infection classes including chronic release/budding
M9 = M2 + ['gam'] #One infection class and one recovery class
M10 = M5 + ['gam'] #Two infected classes and one recovery class
M11 = ['mum']
M12 = M11 +['aff'] #model that accounts for curvature
#mods = [M1, M2, M7, M9, M10]
mods = [M12]
#labs = ['M1', 'M2', 'M7', 'M9', 'M10'] 
labs = ['M11']
labloc = np.arange(len(labs))  

### Set up Results Dataframe ####
allTraits = pd.DataFrame(columns =['OrigRef', 'ModelNumber','mu', 'phi', 'beta', 'lambda','AdjR2', 'AIC'])
allTraits = allTraits[['OrigRef', 'ModelNumber', 'mu', 'phi', 'beta', 'lambda','AdjR2', 'AIC']]


#### Run the models ####

print(' ')
print('###############################')
print('Do fitting')
print('###############################')
for (inf,tag) in zip(exp_set_all,vals_all): #place to select the file we will run MCMC on
    phi_all,beta_all=r_[[]],r_[[]]
    pmod = all_mods(inf,mods[0], nits=1000,pits=100,burnin=500)
    figs,axes = pmod.gen_figs(tag) #crashing here
    allmods = []
    for (mod,lab) in zip(mods,labs):
        print(' ')
        print('dataset: ',tag)
        print('model label: ', lab)
        print('model params: ', mod)
        print(' ')
        pmod = all_mods(inf,mod, nits=1000,pits=100,burnin=500)
        pmod.modellabel = mod
        pmod.do_fitting(axes[0])
        print('AIC: ', pmod.AIC)
        print(' ')
        print('###############################')
        axes[1].plot(pmod.iterations,pmod.likelihoods)
        allmods.append(pmod)
        # phi_all = append(phi_all,exp(pmod.pall[1]))
        allTraits = allTraits.append({'OrigRef': tag, 'mu': pmod.pms['mum'], 'AdjR2': pmod.adj_rsquared, 'AIC': pmod.AIC, "ModelNumber": lab}, ignore_index= True)
        # for a in axes[0]:
        #     a.semilogy()
    axes[0].semilogy()
    mums = r_[[pmod.pms['mum'] for pmod in allmods]]
    # phis = r_[[pmod.pms['phi'] for pmod in allmods]]
    # bets = r_[[pmod.pms['beta'] for pmod in allmods]]
    # dels = r_[[pmod.pms['lambd'] for pmod in allmods]]
    muss = r_[[pmod.pss['mum'] for pmod in allmods]]
    # phss = r_[[pmod.pss['phi'] for pmod in allmods]]
    # bess = r_[[pmod.pss['beta'] for pmod in allmods]]
    # dess = r_[[pmod.pss['lambd'] for pmod in allmods]]
    means,stds = [mums],[muss]
    axes[2][0].scatter(arange(len(mods)),r_[[mod.adj_rsquared for mod in allmods]])
    axes[2][1].scatter(arange(len(mods)),r_[[mod.AIC for mod in allmods]])
    for a in axes[2]:
        a.set_xticks(labloc)
        a.set_xticklabels(labs, fontweight = "bold")
        a.tick_params(direction = "in")
    for (m,e,a) in zip(means,stds,axes[3]):
        a.errorbar(range(m.shape[0]),m,yerr=e,fmt='o')
        a.set_xticks(labloc)
        a.set_xticklabels(labs, fontweight = "bold")
        a.tick_params(direction = "in")
    for f in figs:
        pdf.savefig(f)
        close(f)
pdf.close()

##### Finding the best fits for each study #####

bestTraits = pd.DataFrame(columns =['OrigRef', 'ModelNumber','mu', 'phi', 'beta', 'lambda','AdjR2', 'AIC'])
bestTraits = bestTraits[['OrigRef', 'ModelNumber', 'mu', 'phi', 'beta', 'lambda','AdjR2', 'AIC']]
for tag in file_types:
    ModStudy = allTraits.loc[allTraits.OrigRef == tag]
    if any(any(np.isfinite(ModStudy['AIC']))==True):
        bestModel = ModStudy.loc[ModStudy.AIC == ModStudy.AIC.min()]
        bestTraits = bestTraits.append(bestModel)
    else: 
        bestModel = ModStudy.loc[ModStudy.AdjR2 == ModStudy.AdjR2.max()]
        bestTraits = bestTraits.append(bestModel)
         
writer = pd.ExcelWriter('AlgaeModResults.xlsx', engine='xlsxwriter')

bestTraits.to_excel(writer, sheet_name = "BestFits", index = False)
allTraits.to_excel(writer, sheet_name = "AllFits", index = False)

writer.save()
 
