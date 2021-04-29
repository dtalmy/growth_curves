import pylab as py
import pandas as pd
import scipy
import argparse

def main(VALver):
    #########################################################
    ## Single model run based on tids value
    #########################################################
    # read data
    print(VALver)
    
    master_df = pd.read_csv('../data/input/processed/processed_data.csv',index_col='id')
    abiotic_treatment_df = pd.read_csv('../data/input/preprocessed/reu_2019/treatments.csv',index_col='id')
    abiotic_treatment_df = abiotic_treatment_df[abiotic_treatment_df['treatment']=='Replete']
    master_df = pd.concat((master_df,abiotic_treatment_df))
    treatments = master_df.query('control==False').copy() # remove controls
    tids = treatments.index.unique() # unique ids
    dest_file= 'Gen_sub'+str(VALver)+'.sh'
    print(dest_file)
    with open(dest_file, 'w') as writer:
        writer.write('#!/bin/bash\n')
        writer.write('declare -a arr=(\n')
        writer.write('\n')
        writer.write('')
        writer.write('for i in \\"${arr[@]}\\"\n')
        writer.write('do\n')
        writer.write('  echo \\"$i\\"\n')
        writer.write('   qsub qsub_tids.sh  -v VALtids=$i,VALoutpath=/lustre/haven/proj/UTK0105/Python_runs/aaa  \n')
        writer.write('done\n')
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='generates sub script')
    parser.add_argument('ver', type=str, help='ver tag')
    
    args = parser.parse_args()
    print('ver')
    print(args.ver)
    main(str(args.ver))
