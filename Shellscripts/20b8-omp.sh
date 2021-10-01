#!/bin/bash
#SBATCH --job-name=20b8omp     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p skylake
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --time=01:15:00              # Time limit hrs:min:sec
#SBATCH --output=20b8-omp%j.log     # Standard output and error log
pwd; hostname; date

echo "Running dataset20b8B large dataset on $SLURM_CPUS_ON_NODE CPU cores"
echo "RES _______ 4 __________ THR _________ 8,16,32 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp8x4 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp16x4 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp32x4 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
echo "RES _______ 8 __________ THR _________ 4,8, 16, 32 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp4x8 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp8x8 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp16x8 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp32x8 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
echo "RES _______ 16 __________ THR _________ 2, 4,8, 16 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp2x16 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp4x16 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp8x16 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_omp16x16 -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025

date
