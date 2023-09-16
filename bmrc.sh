#!/bin/bash

#SBATCH -A moru-batty.prj
#SBATCH -D /well/moru-batty/users/gka292/Primaquine-Challenge
#SBATCH -J res
#SBATCH -n 4
#SBATCH -o /well/moru-batty/users/gka292/Primaquine-Challenge/o_and_e_files/output.o%A_%a.out
#SBATCH -e /well/moru-batty/users/gka292/Primaquine-Challenge/o_and_e_files/output.e%A_%a.out
#SBATCH -p long
#SBATCH --array 1-41


echo started=`date`
module purge
module load R/4.1.2-foss-2021b


echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`
Rscript /well/moru-batty/users/gka292/Primaquine-Challenge/run_pop_model.R ${SLURM_ARRAY_TASK_ID} --no-save --no-restore
echo "finished="`date`
exit 0