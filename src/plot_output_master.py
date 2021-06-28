import pandas as pd
import seaborn as sns
import numpy as np
from helpers import *
from collections import Counter

########################################################
# read data
########################################################

# first load subset of data that contain replicate
ndf = pd.read_csv('../data/input/preprocessed/nissimov/with_reps.csv').set_index('organism')
ndf = ndf[ndf['virusName']=='EhV-99B1']
ndfh,ndfv = ndf.loc['H'],ndf.loc['V']

# load master experimental dataset
master_df = pd.read_csv('../data/input/processed/processed_data.csv',index_col='id')
abiotic_treatment_df = pd.read_csv('../data/input/preprocessed/reu_2019/treatments.csv',index_col='id')
abiotic_treatment_df = abiotic_treatment_df[abiotic_treatment_df['treatment']=='Replete']
master_df = pd.concat((master_df,abiotic_treatment_df))
treatments = master_df.query('control==False').copy() # remove controls
tids = treatments.index.unique() # unique ids

# sort experimental data into a single dictionary
datasets = {}
for tid in tids[:-1]:
    df = treatments.loc[tid].copy()
    df.loc[:,'log_sigma'] = 0.2
    df.loc[df.organism == 'H', 'time'] = df.loc[df.organism == 'H', 'time'].copy() -\
                min(df.loc[df.organism == 'H', 'time'])
    df.loc[df.organism == 'V', 'time'] = df.loc[df.organism == 'V', 'time'].copy() -\
                min(df.loc[df.organism == 'V', 'time'])
    datasets[tid] = df

# load posteriors into a dictionary
posteriors = {}
i = 1
print('load data')
for tid in tids[:-1]:
    print(i,tid)
    i = i+1
    #f = '../data/output/jun25th2021/'+tid+'.csv'
    f = '../data/output/jun27th2021_million/'+tid+'.csv'
    n = 100
    num_lines = sum(1 for l in open(f))
    skip_idx = [x for x in range(1, num_lines) if x % n != 0]
    posteriors[tid] = pd.read_csv(f,low_memory=False,skiprows=skip_idx)

########################################################
# load priors
########################################################

mu_prior,phi_prior,beta_prior,tau_prior,H0_prior,V0_prior = load_priors(df)

########################################################
# setup figures
########################################################

f1,ax1 = py.subplots(1,2,figsize=[9,4])
f2,ax2 = py.subplots(3,1,figsize=[12,12])
f3a,ax3a = py.subplots(4,4,figsize=[18,18])
f3b,ax3b = py.subplots(4,4,figsize=[18,18])
f3c,ax3c = py.subplots(4,4,figsize=[18,18])
f3d,ax3d = py.subplots(4,4,figsize=[18,18])
f4,ax4 = py.subplots()

f3a.subplots_adjust(wspace=0.4,hspace=0.6)
f3b.subplots_adjust(wspace=0.4,hspace=0.6)
f3c.subplots_adjust(wspace=0.4,hspace=0.6)
f3d.subplots_adjust(wspace=0.4,hspace=0.6)

# figure font size
fs = 12

########################################################
# do plotting
########################################################

# figure 1
ax1[0].plot(ndfh.time,ndfh.rep1,'-o',label='rep 1')
ax1[0].plot(ndfh.time,ndfh.rep2,'-^',label='rep 2')
ax1[0].plot(ndfh.time,ndfh.rep3,'-*',label='rep 3')
ax1[1].plot(ndfv.time,ndfv.rep1,'-o')
ax1[1].plot(ndfv.time,ndfv.rep2,'-^')
ax1[1].plot(ndfv.time,ndfv.rep3,'-*')
l = ax1[0].legend(prop={'size':fs})
l.draw_frame(False)
ax1[0].semilogy()
ax1[1].semilogy()
for a in ax1:
        a.set_xlabel('Time (hours)',fontsize=fs)
ax1[0].set_ylabel('Host (ml$^{-1}$)',fontsize=fs)
ax1[1].set_ylabel('Virus (ml$^{-1}$)',fontsize=fs)
f1.subplots_adjust(wspace=0.3)

# figure 3
ax3all = np.concatenate((ax3a.flatten(),ax3b.flatten(),ax3c.flatten(),ax3d .flatten()))
hosts = ax3all[0:-1:2]
viruses = ax3all[1:-1:2]
i = 1
allbests = pd.DataFrame()
bestmodels = []
print('begin plotting')
for (hax,vax,tid) in zip(hosts,viruses,tids[:-1]):
    print(i,tid)
    i = i+1
    ddf = datasets[tid]
    models = get_models(ddf)
    hdat =ddf[ddf.organism=='H']
    vdat =ddf[ddf.organism=='V']
    hax.errorbar(hdat.time,np.log(hdat.abundance),yerr=hdat.log_sigma)
    vax.errorbar(vdat.time,np.log(vdat.abundance),yerr=vdat.log_sigma)
    mdf = posteriors[tid] # posteriors
    mi = mdf.loc[mdf.chi==min(mdf.chi)].index[0]
    bestmodelstring = mdf.iloc[mi]['Unnamed: 0']
    bestmodels.append(bestmodelstring)
    bestmodel = models[bestmodelstring]
    bestmodelposteriors = mdf[mdf['Unnamed: 0']==bestmodelstring]
    bestmodelposteriors['algalHost'] = hdat.algalHost.unique()[0]
    bestmodelposteriors['algalHostTaxon'] = hdat.algalHostTaxon.unique()[0]
    bestmodelposteriors['virusName'] = vdat.virusName.unique()[0]
    allbests = pd.concat((allbests,bestmodelposteriors))
    set_optimal_parameters(bestmodel,bestmodelposteriors)
    mod = bestmodel.integrate()
    hax.plot(bestmodel.times,np.log(mod['H']),c='r',lw=2,zorder=2)
    vax.plot(bestmodel.times,np.log(mod['V']),c='r',lw=2,zorder=2)
    hax.set_title(hdat.algalHost.unique()[0])
    vax.set_title(vdat.virusName.unique()[0])
    for a in range(1000):
        set_random_param(bestmodel,bestmodelposteriors)
        mod = bestmodel.integrate()
        hax.plot(bestmodel.times,np.log(mod['H']),c=str(0.8),lw=1,zorder=1)
        vax.plot(bestmodel.times,np.log(mod['V']),c=str(0.8),lw=1,zorder=1)

sns.boxplot(x='algalHostTaxon',y='mu',data=allbests,ax=ax2[0])
sns.boxplot(x='algalHostTaxon',y='phi',data=allbests,ax=ax2[1])
sns.boxplot(x='algalHostTaxon',y='beta',data=allbests,ax=ax2[2])

for a in ax2:
    a.semilogy()

keys, counts = np.unique(bestmodels, return_counts=True)
ax4.bar(keys, counts)
ax4.set_xlabel('Number of infection states')
ax4.set_ylabel('Number of datasets')

########################################################
# save figures
########################################################

f1.savefig('../figures/figure1',bbox_inches='tight',pad_inches=0.1)
f2.savefig('../figures/figure2',bbox_inches='tight',pad_inches=0.1)
f3a.savefig('../figures/figure3a',bbox_inches='tight',pad_inches=0.1)
f3b.savefig('../figures/figure3b',bbox_inches='tight',pad_inches=0.1)
f3c.savefig('../figures/figure3c',bbox_inches='tight',pad_inches=0.1)
f3d.savefig('../figures/figure3d',bbox_inches='tight',pad_inches=0.1)
f4.savefig('../figures/figure4',bbox_inches='tight',pad_inches=0.1)
