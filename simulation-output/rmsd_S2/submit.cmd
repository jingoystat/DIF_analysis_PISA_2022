#!/bin/bash
#####################################################################################
###                                                                                 #
### slurm-mpi.cmd :					                            #				 
### A SLURM submission script for running MPI program in HPC2021 system             #
###                                                                                 #
### Usage:                                                                          #
###    cd to directory containing MPI program executable, then:                     #
###    sbatch <location of this script>/slurm-mpi.cmd                               #    
###                                                                                 #
###    You may modify the script to select different MPI libraries                  #
###      module load impi                 ( Intel MPI Libraries)                    #
###      module load openmpi              ( OpenMPI Libraries)                      #
###                                                                                 #
### - Written by Lilian Chan, HKU ITS (2021-3-2)                                    #
###                                                                                 #
#####################################################################################

#SBATCH --job-name=t-mpi                          # Job name
##SBATCH --mail-user=jingoy@hku.hk                   # Update your email address   
##SBATCH --partition=intel                        # Specific Partition (intel/amd)
##SBATCH --qos=normal                             # Specific QOS (debug/normal/long)
#SBATCH --time=24:00:00                              # Wall time limit (days-hrs:min:sec)
#SBATCH --nodes=1                                 # Total number of compute node(s)
#SBATCH --ntasks=1                               # Total number of MPI tasks (i.e. processes)
#SBATCH --ntasks-per-node=1                      # Number of MPI tasks on each node
#SBATCH --output=%x.out.%j                        # Standard output file  
#SBATCH --error=%x.err.%j                         # Standard error file  

#####################################################################################
### The following stuff will be executed in the first allocated node.               #
### Please don't modify it                                                          #
#####################################################################################
echo "SLURM_NTASKS          : $SLURM_NTASKS"
echo "SLURM_JOB_NUM_NODES   : $SLURM_JOB_NUM_NODES"
echo "SLURM_CPUS_PER_TASK   : $SLURM_CPUS_PER_TASK"
echo "SLURM_CPUS_ON_NODE    : $SLURM_CPUS_ON_NODE"
###################################################################################
cd ${SLURM_SUBMIT_DIR}
OUTFILE=${SLURM_JOB_NAME}.${SLURM_JOBID}
echo ===========================================================
echo "Job Start Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

module load R/4.1.2
R CMD BATCH "--args $arg1 $arg2 $arg3 $arg4 $arg5 $arg6" /home/jingoy/dif_lrt_J15_p10_random/main.R


echo "Job Finish Time is `date "+%Y/%m/%d -- %H:%M:%S"`"

exit 0
