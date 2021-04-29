import pylab as py
import pandas as pd
import scipy
import argparse
import subprocess
import os

def main(VALver):
    #########################################################
    ## Single model run based on tids value
    #########################################################
    # read data
    print(VALver)
    dest_file= 'Gen_sub'+str(VALver)+'.sh'
    print(dest_file)
    #
    master_df = pd.read_csv('../data/input/processed/processed_data.csv',index_col='id')
    abiotic_treatment_df = pd.read_csv('../data/input/preprocessed/reu_2019/treatments.csv',index_col='id')
    abiotic_treatment_df = abiotic_treatment_df[abiotic_treatment_df['treatment']=='Replete']
    master_df = pd.concat((master_df,abiotic_treatment_df))
    treatments = master_df.query('control==False').copy() # remove controls
    tids = treatments.index.unique() # unique ids
    
    with open(dest_file, 'w') as writer:
        writer.write('#!/bin/bash\n')
        writer.write('declare -a arr=(\n')
        for did in tids:
            print(did)
            writer.write('\"'+did+'\"\n')
        writer.write(')\n')
        writer.write('for i in \"${arr[@]}\"\n')
        writer.write('do\n')
        writer.write('   echo \"$i\"\n')
        writer.write('   qsub qsub_tids.sh  -v VALtids=$i,VALoutpath=/lustre/haven/proj/UTK0105/Python_runs/aaa  \n')
        writer.write('done')
    result=subprocess.run(["ls", "-l"], capture_output=True, text=True)    
    print(result.stdout)
    cmdresult=subprocess.run(["touch "+dest_file], capture_output=True, text=True)
    print(cmdresult.stdout)
    print(cmdresult.stderr)
    
    os.chmod(dest_file, stat.S_IRUSR | stat.S_IWUSR |stat.S_IXUSR |stat.S_IRGRP | stat.S_IROTH)
    
    result=subprocess.run(["ls", "-l"], capture_output=True, text=True)    
    print(result.stdout)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='generates sub script')
    parser.add_argument('ver', type=str, help='ver tag')
    
    args = parser.parse_args()
    print('ver')
    print(args.ver)
    main(str(args.ver))
