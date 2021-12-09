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
    did = VALtids
    print(did)
    main_df = get_master_dataframe()
    df = main_df.loc[did].copy()
    tpdf = fit_all(df) # here is where the main work is done
    tpdf.close()
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='run model for TIDS')
    parser.add_argument('tids', type=str, help='the tids name')
    parser.add_argument('outpath', type=str, help='path to save output')
    args = parser.parse_args()
    print('tids, outpath')
    print(args.tids, args.outpath, sep=',')
    main(str(args.tids), str(args.outpath))
    
    
