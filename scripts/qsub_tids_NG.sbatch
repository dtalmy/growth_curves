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
OD_base='/lustre/isaac/proj/UTK0105/GrowthCurves/Test/growth_curves'
JOBID=$SLURM_JOB_ID
echo JOBID : $JOBID
echo VALtids   : $VALtids
echo VALoutpath   : $VALoutpath
#Setup run 
cd $OD_base
mkdir $OD_base/$JOBID
cd $OD_base/$JOBID
git clone -b master https://github.com/dtalmy/growth_curves
cd $OD_base/$JOBID/growth_curves/src 
#
#run mo
# mail/move/delete
#mailx -a '../figures/'$PBS_JOBID'_test_selection.pdf' -s 'RESULTS_'$VALpdfname ecarr@utk.edu < /dev/null
cp ../data/output/* $OD_base/aaa/posteriors/
cp ../figures/* $OD_base/aaa/figures/
cp ../data/params/final/* $OD_base/aaa/params/
#make DIR readable by group
chmod +R g+r $OD_base/aaa
chmod +R g+r $OD_base/$JOBID
chmod +R g+r $OD_base/job_data
#cd /lustre/haven/proj/UTK0105/Python_runs
#rm -rf ./$PBS_JOBID
