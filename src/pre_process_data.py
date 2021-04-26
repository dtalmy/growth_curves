'''
Data Formatting for Host-Virus Interactions Data

Created @audralhinson
Created Aug 16th, 2019
Updated Sept 10th, 2019
'''

from pylab import *
import pandas as pd
import glob

# function to process data
def read_txt_file(fileloc):
    filename = re.sub(path,'',fileloc)
    tag = filename.split('_')[0]
    print('reading ',filename)
    if 'StdDev' not in str(filename): 
        HVDat = pd.read_csv(fileloc,names = ['time','abundance'])
    else: 
        HVDat = pd.read_csv(fileloc, names = ['time','abundance', 'lowlim','highlim'])
        if tag not in ('Nagasaki29','Goa41','Ou55','Tomaru39'): # some of the datasets report inappropriate error
            HVDat['stddev']= (HVDat['lowlim'] + HVDat['highlim']) / 2.0 # if error are log-transformed, keep them that way
            HVDat['stddev'] = HVDat['stddev'].apply(lambda x: round(x,2))
        HVDat = HVDat.drop(['lowlim','highlim'],axis=1)
    if 'Time' in str(filename):
        HVDat['time'] = HVDat.apply(lambda x: x*24)
    HVDat['time'] = HVDat['time'].apply(lambda x: round(float(x),2))
    if 'Log' in str(filename):
        HVDat['abundance'] = pow(10, HVDat['abundance'])
        HVDat['logscale'] = True
    else:
        HVDat['logscale'] = False
    if ('Control' in str(filename)) or ('Infected' in str(filename)):
        HVDat['organism'] = 'H'
        if ('Control' in str(filename)):
            HVDat['control'] = True
        else:
            HVDat['control'] = False
    else:
        HVDat['organism'] = 'V'
        HVDat['control'] = False
    HVDat['id'] = tag
    return (HVDat)

# specify directories and files where data are
path = '../data/input/preprocessed/audra/'
sourceInfo = pd.read_csv(path+'SourceTableInfo.csv', header = 0)
Allfiles = glob.glob(path+'*.txt')

# read data into dataframe
raw_df = pd.DataFrame()
for ind,tag in enumerate(Allfiles):
    temp_df = read_txt_file(tag)
    raw_df = raw_df.append(temp_df)

# merge raw data with medadata file
master_df = pd.merge(raw_df, sourceInfo, on ='id')
master_df = master_df.set_index('id')
master_df.to_csv('../data/input/processed/'+'processed_data.csv')
