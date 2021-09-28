#!/bin/bash
#SBATCH --job-name=17gbhyb      # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p skylake
#SBATCH --nodes=2                    # Run all processes on a single node	
#SBATCH --time=01:15:00              # Time limit hrs:min:sec
#SBATCH --output=17gb-hyb%j.log     # Standard output and error log
pwd; hostname; date

echo "Running 17GB large dataset on $SLURM_CPUS_ON_NODE CPU cores"
mpirun -n 2 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_hyb -f /mnt/beegfs/smdupor/1.7g.txt -n 4 -r 1659025

date
