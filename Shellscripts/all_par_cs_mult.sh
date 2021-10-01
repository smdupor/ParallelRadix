#!/bin/bash
#SBATCH --job-name=frParSer     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=smdupor@ncsu.edu    # Where to send mail	
#SBATCH -p skylake
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --time=12:15:00              # Time limit hrs:min:sec
#SBATCH --output=frParSer-%j.log     # Standard output and error log
pwd; hostname; date

echo "Running PARALLEL COUNT Sort for All dataset on $SLURM_CPUS_ON_NODE CPU cores"

echo "64m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/64m.txt -n 2 -r 1659025

echo "84m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/84m.txt -n 2 -r 1659025

echo "479m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/479m.txt -n 2 -r 1659025

echo "1.7 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/1.7g.txt -n 2 -r 1659025

echo "16b8m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset16b8.txt -n 2 -r 1659025

echo "20b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset20b8.txt -n 2 -r 1659025

echo "24b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset24b8.txt -n 2 -r 1659025

echo "28b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset28b8.txt -n 2 -r 1659025

echo "30b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset30b8.txt -n 2 -r 1659025

echo "13gb ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/13gb.txt -n 2 -r 1659025

echo "-------------------------------------------------- 4 CORES-----------------"
echo "64m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/64m.txt -n 4 -r 1659025

echo "84m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/84m.txt -n 4 -r 1659025

echo "479m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/479m.txt -n 4 -r 1659025

echo "1.7 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/1.7g.txt -n 4 -r 1659025

echo "16b8m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset16b8.txt -n 4 -r 1659025

echo "20b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset20b8.txt -n 4 -r 1659025

echo "24b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset24b8.txt -n 4 -r 1659025

echo "28b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset28b8.txt -n 4 -r 1659025

echo "30b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset30b8.txt -n 4 -r 1659025

echo "13gb ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/13gb.txt -n 4 -r 1659025

echo "-------------------------------------------------- 8 CORES-----------------"

echo "64m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/64m.txt -n 8 -r 1659025

echo "84m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/84m.txt -n 8 -r 1659025

echo "479m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/479m.txt -n 8 -r 1659025

echo "1.7 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/1.7g.txt -n 8 -r 1659025

echo "16b8m ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset16b8.txt -n 8 -r 1659025

echo "20b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset20b8.txt -n 8 -r 1659025

echo "24b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset24b8.txt -n 8 -r 1659025

echo "28b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset28b8.txt -n 8 -r 1659025

echo "30b8 ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/dataset30b8.txt -n 8 -r 1659025

echo "13gb ___________"
mpirun -n 1 -bootstrap slurm /home/smdupor/ParallelRadix/bin/main_ser_rad -f /mnt/beegfs/smdupor/13gb.txt -n 8 -r 1659025

date
