import pylab as py
import pandas as pd
import scipy
from helpers import *
from matplotlib.backends.backend_pdf import PdfPages
import ODElib

# read data
master_df = pd.read_csv('../data/input/processed/processed_data.csv',index_col='id')
abiotic_treatment_df = pd.read_csv('../data/input/preprocessed/reu_2019/treatments.csv',index_col='id')
abiotic_treatment_df = abiotic_treatment_df[abiotic_treatment_df['treatment']=='Replete']
master_df = pd.concat((master_df,abiotic_treatment_df))
treatments = master_df.query('control==False').copy() # remove controls
tids = treatments.index.unique() # unique ids

# loop over each datset and save model output to pdf
for did in tids:
    print(did)
    df = treatments.loc[did].copy()
    df.loc[:,'log_sigma'] = 0.2 # define uncertainty in data 
    df.loc[df.organism == 'H', 'time'] = df.loc[df.organism == 'H', 'time'].copy() -\
                min(df.loc[df.organism == 'H', 'time'])
    df.loc[df.organism == 'V', 'time'] = df.loc[df.organism == 'V', 'time'].copy() -\
                min(df.loc[df.organism == 'V', 'time'])
    tpdf = fit_all(df)
    tpdf.close()
