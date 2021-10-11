#!/bin/bash
#SBATCH --job-name=30b8hyb      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p normal
#SBATCH --nodes=2                    # Run all processes on a single node	
#SBATCH --time=01:15:00              # Time limit hrs:min:sec
#SBATCH --output=30b8-hyb%j.log     # Standard output and error log
pwd; hostname; date

echo "Running 30b8GB large dataset on $SLURM_CPUS_ON_NODE CPU cores"
mpirun -n 2 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_hyb -f /mnt/beegfs/smdupor/dataset30b8.txt -n 4 -r 1659025

date


