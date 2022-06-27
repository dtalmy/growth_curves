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
import pathlib
import time  


def main(VALver):
    #########################################################
    ## Single model run based on tids value
    #########################################################
    # Create file name
    #
    dest_file= 'Sub_TIDS'+str(VALver)+'_NG.sh'
    PD_cwd=pathlib.Path.cwd()
    print(dest_file)
    print('\n')
    ####
    # OD_base is where the individual TIDS will run and the combined data will be saved.
    # edit basepath to control where the run will go
    #
    # setup a unique time tag, it does not need to be epoch time but gurantees unique and timestamped
    # on linux terminal date -d @$epochtime convert back to date
    # converting to int drop microseconds without floor
    epochtime= str(int(time.time()))
    basepath=str('/lustre/isaac/proj/UTK0105/GrowthCurves/growth_curves_')+epochtime
    #
    OD_base=pathlib.Path(basepath)
    OD_base.mkdir(parents=True, exist_ok=True)
    # The helper functions have hard coded relative paths. It expects to be in scripts
    OD_subout=str('combined')
    pl_combined=OD_base / OD_subout
    pl_combined.mkdir(parents=True, exist_ok=True)
    OD_figures=str('figures')
    pl_figures=pl_combined / OD_figures
    pl_figures.mkdir(parents=True, exist_ok=True)
    OD_params=str('params')
    pl_params=pl_combined / OD_params
    pl_params.mkdir(parents=True, exist_ok=True)
    #posteriors
    OD_posteriors=str('posteriors')
    pl_posteriors=pl_combined / OD_posteriors
    pl_posteriors.mkdir(parents=True, exist_ok=True)
    #
    #
    #filename=str('Dsum_LineSum.png')
    
    OF_dest=PD_cwd /str(dest_file)
    print('Out dir set to: ' + OD_base.as_posix(), file = sys.stdout )
    print('Out file script: ' + OF_dest.as_posix(), file = sys.stdout )
    
    
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
    with open(OF_dest.as_posix(), 'w') as writer:
        writer.write('#!/bin/bash\n')
        writer.write('declare -a arr=(\n')
        for did in tids:
            writer.write('\"'+did+'\"\n')
        writer.write(')\n')
        writer.write('for i in \"${arr[@]}\"\n')
        writer.write('do\n')
        writer.write('   echo \"$i\"\n')
        writer.write(str('   sbatch --export=VALtids=$i,OD_base=\"')+OD_base.as_posix()+str('\"  qsub_tids_NG.sbatch  \n') )
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
