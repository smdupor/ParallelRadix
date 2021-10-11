#!/bin/bash
#SBATCH --job-name=20b8mp     # Job name
#SBATCH --mail-type=FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p normal
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --time=10:15:00              # Time limit hrs:min:sec
#SBATCH --output=20b8-omp%j.log     # Standard output and error log
pwd; hostname; date

echo "Running dataset20b8B large dataset on $SLURM_CPUS_ON_NODE CPU cores"
echo "RES _______ 4 __________________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParSub/bin/run-graph-openmp -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1
echo "RES _______ 4 __________________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParSub/bin/run-graph-openmp -f /mnt/beegfs/smdupor/dataset20b8.txt -n 8 -r 1
echo "RES _______ 4 __________________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParSub/bin/run-graph-openmp -f /mnt/beegfs/smdupor/dataset20b8.txt -n 16 -r 1

date
