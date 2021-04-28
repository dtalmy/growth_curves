import pylab as py
import pandas as pd
import scipy
from helpers import *
from matplotlib.backends.backend_pdf import PdfPages
import ODElib
import argparse

def main(VALtids,DATAdir):
    #########################################################
    ## Single model run based on tids value
    #########################################################
    # read data
    print(VALtids)
    print(DATAdir)
    master_df = pd.read_csv('../data/input/processed/processed_data.csv',index_col='id')
    abiotic_treatment_df = pd.read_csv('../data/input/preprocessed/reu_2019/treatments.csv',index_col='id')
    abiotic_treatment_df = abiotic_treatment_df[abiotic_treatment_df['treatment']=='Replete']
    master_df = pd.concat((master_df,abiotic_treatment_df))
    treatments = master_df.query('control==False').copy() # remove controls
    #tids = treatments.index.unique() # unique ids

    # process tids value in datset and save model output to pdf
    did = VALtids
    print(did)
    df = treatments.loc[did].copy()
    df.loc[:,'log_sigma'] = 0.2 # define uncertainty in data 
    df.loc[df.organism == 'H', 'time'] = df.loc[df.organism == 'H', 'time'].copy() -\
                min(df.loc[df.organism == 'H', 'time']) # remove non-zero initial timepoints
    df.loc[df.organism == 'V', 'time'] = df.loc[df.organism == 'V', 'time'].copy() -\
                min(df.loc[df.organism == 'V', 'time']) # same for virus
    tpdf = fit_all(df) # here is where the main work is done
    tpdf.close()
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='run model for TIDS')
    parser.add_argument('--tids', metavar='path', required=True,
                        help='the tids name')
    parser.add_argument('--outpath', metavar='path', required=True,
                        help='path to save output')
    args = parser.parse_args()
    print('tids, path')
    print(args.tids, args.path, sep=',')
    main(str(args.tids), str(args.outpath))
    
    
