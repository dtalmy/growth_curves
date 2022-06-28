#!/bin/bash -x
#SBATCH -J TIDSsub
#SBATCH -A ACF-UTK0105
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=/lustre/isaac/proj/UTK0105/logfiles/GrowCurve.%j.out
#SBATCH --error=/lustre/isaac/proj/UTK0105/logfiles/GrowCurve.%j.err
#SBATCH --qos=condo
#SBATCH --partition=condo-dtalmy
echo ---Modules 
module purge
#module load git
#load anaconda
module load anaconda3/
source $ANACONDA3_SH
conda activate PY39
python3 -V
module list
# od base is passed on the submission line 
#OD_base='/lustre/isaac/proj/UTK0105/GrowthCurves/Test/growth_curves'
JOBID=$SLURM_JOB_ID
echo JOBID : $JOBID
echo VALtids   : $VALtids
echo OD_base   : $OD_base
#
# chk dir and make if needed.  This is overly redundent but if you edit submission script you would need to create these. Also created in the generation script.
####
# OD_base is where the individual TIDS will run and the combined data will be saved.
#
#OD_base=pathlib.Path('/lustre/isaac/proj/UTK0105/GrowthCurves/growth_curves/')
# The helper functions have hard coded relative paths. It expects to be in scripts
mkdir -p $OD_base
OD_subout='combined'
pl_combined=$OD_base'/'$OD_subout
mkdir -p $pl_combined
OD_figures='figures'
pl_figures=$pl_combined'/'$OD_figures
mkdir -p $pl_figures
OD_params='params'
pl_params=$pl_combined'/'$OD_params
mkdir -p $pl_params
#posteriors
OD_posteriors='posteriors'
pl_posteriors=$pl_combined'/'$OD_posteriors
mkdir -p $pl_posteriors
#
#######################
#Setup run 
# jobid will be unique for each tids run in the set 
cd $OD_base
mkdir $OD_base'/'$JOBID
cd $OD_base/$JOBID
git clone -b master https://github.com/dtalmy/growth_curves
###
# if for some reson you want to copy a dir instead of get it from github
# Model_Dir='/lustre/isaac/proj/UTK0105/GrowthCurves/Test/growth_curves'
# cp -r $$Model_Dir $OD_base/$JOBID
cd $OD_base/$JOBID/growth_curves/src 
#
#run model
VALpdfname="tids_"$VALtids
echo VALpdfname : $VALpdfname
python3 single_tids_model.py $VALtids '../figures/'
# mail/move/delete
#mailx -a '../figures/'$PBS_JOBID'_test_selection.pdf' -s 'RESULTS_'$VALpdfname ecarr@utk.edu < /dev/null
#
# copy output
#
cp ../data/output/* $OD_base/combined/posteriors/
cp ../figures/* $OD_base/combined/figures/
cp ../data/params/final/* $OD_base/combined/params/
#make DIR readable by group
chmod +R g+r $OD_base/aaa
chmod +R g+r $OD_base/$JOBID
chmod +R g+r $OD_base/job_data
#cd /lustre/haven/proj/UTK0105/Python_runs
#rm -rf ./$PBS_JOBID
