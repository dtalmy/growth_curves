import pylab as py
import pandas as pd
import numpy as np
import scipy
import sys
sys.path.append('../src')
from helpers import *
import argparse
import subprocess
import os
import stat

def main(VALver):
    #########################################################
    ## Single model run based on tids value
    #########################################################
    # Create file name
    #
    dest_file= 'Sub_TIDS'+str(VALver)+'.sh'
    print(dest_file)
    print('\n')
    #
    # read data
    #
    main_df = get_master_dataframe()
    tids = main_df.index.unique() # unique ids
    print('Found tids:\n')
    for did in tids:
            print(did)
    #
    # Write bash script
    # the with automatically closes the file
    #
    with open(dest_file, 'w') as writer:
        writer.write('#!/bin/bash\n')
        writer.write('declare -a arr=(\n')
        for did in tids:
            writer.write('\"'+did+'\"\n')
        writer.write(')\n')
        writer.write('for i in \"${arr[@]}\"\n')
        writer.write('do\n')
        writer.write('   echo \"$i\"\n')
        writer.write('   qsub qsub_tids.sh  -v VALtids=$i,VALoutpath=/lustre/haven/proj/UTK0105/Python_runs/aaa  \n')
        writer.write('done')
    #
    # Set file to be executible for both user and group
    #    
    print('\n')
    print("Set file permission for execute")
    result=subprocess.run(["ls", "-l",dest_file], capture_output=True, text=True)    
    print(result.stdout)
        
    os.chmod(dest_file, stat.S_IRUSR | stat.S_IWUSR |stat.S_IXUSR |stat.S_IRGRP |stat.S_IWGRP |stat.S_IXGRP | stat.S_IROTH)
    
    result=subprocess.run(["ls", "-l",dest_file], capture_output=True, text=True)    
    print(result.stdout)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='generates sub script')
    parser.add_argument('ver', type=str, help='ver tag')
    
    args = parser.parse_args()
    print('ver')
    print(args.ver)
    main(str(args.ver))
