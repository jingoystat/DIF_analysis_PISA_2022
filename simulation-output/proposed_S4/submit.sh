#!/bin/bash
   for s in 82 83 84 85 86 87 88 89
    do
    for r in 2
do
for n in 10000 20000
do
for J in 30
do
for p in 20
do

    sbatch --export arg1=$s,arg2=$r,arg3=$n,arg4=$J,arg5=$p -o /dev/null -e /dev/null submit.cmd

done
done
done
done
done


