#!/bin/bash
sbatch ./Shellscripts/64m-omp.sh
sbatch ./Shellscripts/64m-mpi.sh
sbatch ./Shellscripts/64m-hyb.sh
sbatch ./Shellscripts/84m-omp.sh
sbatch ./Shellscripts/84m-mpi.sh
sbatch ./Shellscripts/84m-hyb.sh
sbatch ./Shellscripts/479m-omp.sh
sbatch ./Shellscripts/479m-mpi.sh
sbatch ./Shellscripts/479m-hyb.sh
sbatch ./Shellscripts/1.7gb-omp.sh
sbatch ./Shellscripts/1.7gb-mpi.sh
sbatch ./Shellscripts/1.7gb-hyb.sh
sbatch ./Shellscripts/13gb-omp.sh
sbatch ./Shellscripts/13gb-mpi.sh
sbatch ./Shellscripts/13gb-hyb.sh