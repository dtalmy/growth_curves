#!/bin/bash

declare -a arr=(
"Baudoux32"
"Eissler40"
"Kim52"
"Nagasaki29"
"Tomaru51"
"Slagter2016"
"Baudoux33"
"Gao41"
"Kim71"
"Nagasaki42"
"Tomaru56"
"Slagter2016a"
"Baudoux34"
"Johannessen68"
"Kim72"
"Sandaa26"
"Tomaru57"
"Maat2016a"
"Baudoux3"
"Johannessen69"
"Kimura54"
"Tomaru39"
"Toyoda53"
"Maat2016b"
"Castberg24"
"Johannessen70"
"Kimura66"
"Tomaru43"
"Steenhauer2016"
"Kimura67"
)

for i in "${arr[@]}"
do
   echo "$i"
   qsub qsub_tids.sh  -v VALtids=$i,VALoutpath=/lustre/haven/proj/UTK0105/Python_runs/aaa  
done
