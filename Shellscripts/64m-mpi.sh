#!/bin/bash
#SBATCH --job-name=64mmpi      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p skylake
#SBATCH --nodes=2                    # Run all processes on a single node	
#SBATCH --time=01:15:00              # Time limit hrs:min:sec
#SBATCH --output=64m-mpi%j.log     # Standard output and error log
pwd; hostname; date

echo "Running 64m large dataset on $SLURM_CPUS_ON_NODE CPU cores"
mpirun -n 2 -bootstrap slurm /home/smdupor/ParSub/bin/run-graph-mpi -f /mnt/beegfs/smdupor/64m.txt -n 4 -r 1

date
