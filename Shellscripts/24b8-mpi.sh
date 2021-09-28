#!/bin/bash
#SBATCH --job-name=24b8mpi      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p skylake
#SBATCH --nodes=2                    # Run all processes on a single node	
#SBATCH --time=01:15:00              # Time limit hrs:min:sec
#SBATCH --output=24b8-mpi%j.log     # Standard output and error log
pwd; hostname; date

echo "Running 24b8 large dataset on $SLURM_CPUS_ON_NODE CPU cores"
mpirun -n 2 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_mpi -f /mnt/beegfs/smdupor/dataset24b8.txt -n 4 -r 1659025

date
