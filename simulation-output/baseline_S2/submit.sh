#!/bin/bash
   for s in 89
    do
    for r in 2
do
for n in 20000
do
for J in 15
do
for p in 10
do

    sbatch --export arg1=$s,arg2=$r,arg3=$n,arg4=$J,arg5=$p -o /dev/null -e /dev/null submit.cmd

done
done
done
done
done


