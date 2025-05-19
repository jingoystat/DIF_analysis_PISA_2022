#!/bin/bash
    for s in 80 82 84 86 88
    do
    for r in 2
do
for n in 10000 20000
do
for J in 30
do
for p in 20
do
for t in 0.05 0.1 0.15
do
    sbatch --export arg1=$s,arg2=$r,arg3=$n,arg4=$J,arg5=$p,arg6=$t -o /dev/null -e /dev/null submit.cmd

done
done
done
done
done
done

