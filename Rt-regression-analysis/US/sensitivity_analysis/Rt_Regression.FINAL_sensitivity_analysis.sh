#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 0-08:00              # Runtime in D-HH:MM
#SBATCH -p test,shared,xlin,xlin-lab,hsph,canstat-p01	        # Partition to submit to
#SBATCH --mem=30000            	# Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=xihaoli@g.harvard.edu  # Email to which notifications will be sent

module purge
module load gcc/7.1.0-fasrc01
module load intel-mkl/2017.2.174-fasrc01
module load R/3.5.1-fasrc02
export R_LIBS_USER=$HOME/apps/R-3.5.1-MKL:/n/helmod/apps/centos7/Comp/gcc/7.1.0-fasrc01/R_packages/3.5.1-fasrc02

/n/home01/xihaoli/R-3.5.1/bin/R --vanilla --args ${SLURM_ARRAY_TASK_ID} < $1 > "${1}.${SLURM_ARRAY_TASK_ID}.out" 2> "${1}.${SLURM_ARRAY_TASK_ID}.err"

