import pylab as py
import pandas as pd
import scipy
from helpers import *
from matplotlib.backends.backend_pdf import PdfPages
import ODElib

# read data
main_df = get_master_dataframe()
tids = tids = main_df.index.unique() # get unique indexes
df = main_df.loc[tids[0]].copy() # select one index for testing
tpdf = fit_all_dir(df) # here is where the main work is done
tpdf.close() # save pdf
