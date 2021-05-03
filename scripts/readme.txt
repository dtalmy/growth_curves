You need to check out the repository somewhere to gain access to the scripts.

git clone -b master https://github.com/dtalmy/growth_curves

You can update using

git pull


Go to scripts directory

Qsub submissions must be done on a login node. Not in an interactive session.

Python script to generate Submission script.
Generate_Bash_sub.py
Run the file with the same python3.9 used for model. It will take a tag at the end to  make generated file unique.  The script makes the resulting bash script executable for you.  The bash script is similar to the one described below


[werdna@acf-sk055 scripts]$ python3.9 Generate_Bash_sub.py v2
ver
v2
Sub_TIDSv2.sh


Found tids:

Baudoux32
Baudoux33
Baudoux34
Baudoux3
Castberg24
Eissler40
Gao41
Johannessen68
Johannessen69
Johannessen70
Kim52
Kim71
Kim72
Kimura54
Kimura66
Kimura67
Nagasaki29
Nagasaki42
Sandaa26
Tomaru39
Tomaru43
Tomaru51
Tomaru56
Tomaru57
Toyoda53
Steenhauer2016
Slagter2016
Slagter2016a
Maat2016a
Maat2016b


Set file permission for execute
-rw-r--r-- 1 werdna tug2259 522 Apr 30 11:58 Sub_TIDSv2.sh

-rwxrwxr-- 1 werdna tug2259 522 Apr 30 11:58 Sub_TIDSv2.sh

[werdna@acf-sk055 scripts]$



BASH:
This script is designed to automatically submit the “qsub qsub_tids.sh”. It loops over an array, fills in the arguments, and submits.  Each submission is an independent qsub submission.

Change permission:
[werdna@acf-login8 scripts]$ ls -la
total 20
drwxr-sr-x 2 werdna tug2259 4096 Apr 29 10:29 .
drwxr-sr-x 7 werdna tug2259 4096 Apr 29 10:23 ..
-rw-r--r-- 1 werdna tug2259  525 Apr 29 10:23 bash_qsub_generator.sh
-rw-r--r-- 1 werdna tug2259 2001 Apr 29 10:29 qsub_isaac_serial.sh
-rw-r--r-- 1 werdna tug2259 2019 Apr 29 10:23 qsub_tids.sh
[werdna@acf-login8 scripts]$ chmod ug+x bash_qsub_generator.sh
[werdna@acf-login8 scripts]$ ls -la
total 20
drwxr-sr-x 2 werdna tug2259 4096 Apr 29 10:29 .
drwxr-sr-x 7 werdna tug2259 4096 Apr 29 10:23 ..
-rwxr-xr-- 1 werdna tug2259  525 Apr 29 10:23 bash_qsub_generator.sh
-rw-r--r-- 1 werdna tug2259 2001 Apr 29 10:29 qsub_isaac_serial.sh
-rw-r--r-- 1 werdna tug2259 2019 Apr 29 10:23 qsub_tids.sh

Run:
 ./bash_qsub_generator.sh

[werdna@acf-login8 scripts]$ ./bash_qsub_generator.sh
Baudoux32

Notice: Setting vmem=9832MB due to selected partitions and number of cores.
        If you need more memory, please resubmit your job requesting more cores.
        If you need more memory per core, please resubmit requesting node types with higher memory per core.

Notice: Using default node access policy: shared

4773675.apollo-acf
Eissler40
.
.
.
4773703.apollo-acf
Kimura67

Notice: Setting vmem=9832MB due to selected partitions and number of cores.
        If you need more memory, please resubmit your job requesting more cores.
        If you need more memory per core, please resubmit requesting node types with higher memory per core.

Notice: Using default node access policy: shared

4773704.apollo-acf
[werdna@acf-login8 scripts]$




TIDS:
qsub qsub_tids.sh  -v VALtids=Tomaru56,VALoutpath=/lustre/haven/proj/UTK0105/Python_runs/aaa

[werdna@acf-login8 scripts]$ qsub qsub_tids.sh  -v VALtids=Tomaru56,VALoutpath=/lustre/haven/proj/UTK0105/Python_runs/aaa

Notice: Setting vmem=9832MB due to selected partitions and number of cores.
        If you need more memory, please resubmit your job requesting more cores.
        If you need more memory per core, please resubmit requesting node types with higher memory per core.

Notice: Using default node access policy: shared

4773674.apollo-acf
[werdna@acf-login8 scripts]$ qstat
Job ID                    Name             User            Time Use S Queue
------------------------- ---------------- --------------- -------- - -----
4773656.apollo-acf         STDIN            werdna          00:00:03 R batch
4773673.apollo-acf         Python_serial    werdna          00:00:02 C batch
4773674.apollo-acf         TIDS_submission  werdna                 0 Q batch
[werdna@acf-login8 scripts]$

SERIAL:
qsub qsub_isaac_serial.sh


[werdna@acf-login8 scripts]$ qsub qsub_isaac_serial.sh

Notice: Setting vmem=9832MB due to selected partitions and number of cores.
        If you need more memory, please resubmit your job requesting more cores.
        If you need more memory per core, please resubmit requesting node types with higher memory per core.

Notice: Using default node access policy: shared

4773673.apollo-acf
[werdna@acf-login8 scripts]$ qstat
Job ID                    Name             User            Time Use S Queue
------------------------- ---------------- --------------- -------- - -----
4773656.apollo-acf         STDIN            werdna          00:00:02 R batch
4773673.apollo-acf         Python_serial    werdna                 0 Q batch
[werdna@acf-login8 scripts]$



