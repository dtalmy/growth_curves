#!/bin/bash
#PBS -S /bin/bash
#PBS -A ACF-UTK0105
#PBS -N Python_serial
#PBS -m abe
#PBS -M ecarr@utk.edu
#PBS -o /lustre/haven/proj/UTK0105/Python_runs/job_data/$PBS_JOBID-out.txt
#PBS -e /lustre/haven/proj/UTK0105/Multi_Model/job_data/$PBS_JOBID-err.txt
#PBS -l nodes=1:ppn=2
#PBS -l partition=general
#PBS -l feature=skylake
#PBS -l qos=condo
#PBS -l walltime=96:00:00
#PBS - VALtag test
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo "$PBS_NODEFILE"
echo $(more $PBS_NODEFILE)
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo -----------------------------------------------------
echo -----------------------------------------------------
echo ---Modules 
module load git/2.13.0 
module list
echo PBS_JOBID : $PBS_JOBID
#
#Setup run 
cd /lustre/haven/proj/UTK0105/Python_runs
mkdir ./$PBS_JOBID
cd ./$PBS_JOBID
git clone -b master https://github.com/dtalmy/growth_curves
cd ./growth_curves/src 
#
#run model
# do not use pdfname for this version
#VALpdfname=$PBS_JOBID"_"$VALtag
#echo VALpdfname : $VALpdfname
/lustre/haven/proj/UTK0105/usr/bin/python3.9 test_model.py 
# mail/move/delete
#mailx -a '../figures/'$PBS_JOBID'_test_selection.pdf' -s 'RESULTS_'$VALpdfname ecarr@utk.edu < /dev/null
#cp ../figures/$PBS_JOBID* /lustre/haven/proj/UTK0105/Python_runs/figures/
#rm -rf ./$PBS_JOBID
