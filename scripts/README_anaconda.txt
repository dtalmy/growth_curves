Setup anaconda:
Module purge

module load anaconda3/

source $ANACONDA3_SH


#### to create 

conda create -n PY39 python=3.9 numpy scipy pandas matplotlib



#### now load it and add odelib

conda activate PY39

conda config --add channels conda-forge

conda config --set channel_priority strict

pip3 install git+https://github.com/SEpapoulis/ODElib#egg=ODElib


The PY39 env is referenced in the scripts and now should be able to be loaded during the batch jobs.  There is no reason you cannot use this python 3.9 install for other things.
