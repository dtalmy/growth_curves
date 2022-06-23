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

def main(VALver):
    #########################################################
    ## Single model run based on tids value
    #########################################################
    # Create file name
    #
    dest_file= 'Sub_TIDS'+str(VALver)+'_NG.sh'
    print(dest_file)
    print('\n')
    OD_base=pathlib.Path('/lustre/isaac/proj/UTK0105/GrowthCurves/growth_curves/')
    # The helper functions have hard coded relative paths. It expects to be in scripts
    os.chdir(OD_base / 'scripts')
    OD_subout=str('aaa')
    pl_aaa=OD_base / OD_subout
    pl_aaa.mkdir(parents=True, exist_ok=True)
    OD_figures=str('figures')
    pl_figures=pl_aaa / OD_figures
    pl_figures.mkdir(parents=True, exist_ok=True)
    OD_params=str('params')
    pl_params=pl_aaa / OD_params
    pl_params.mkdir(parents=True, exist_ok=True)
    #posteriors
    OD_posteriors=str('posteriors')
    pl_posteriors=pl_aaa / OD_posteriors
    pl_posteriors.mkdir(parents=True, exist_ok=True)
    #
    #
    filename=str('Dsum_LineSum.png')
    OD_unique=OD_base.joinpath(OD_subout)
    OF_dest=OD_base /'scripts'/str(dest_file)
    print('Out dir set to: ' + OD_unique.as_posix(), file = sys.stdout )
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
        writer.write(str('   sbatch --export=VALtids=$i,VALoutpath=\"')+OD_unique.as_posix()+str('\"  qsub_tids_NG.sbatch  \n') )
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
