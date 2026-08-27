#!/bin/bash
#SBATCH --account=m72-gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=01:00:00
#SBATCH --job-name=w_diag3
#SBATCH --output=%x-%j.out
module load gromacs-amd-gfx90a/2026.1-mixed
export OMP_NUM_THREADS=8
cd /scratch/m72/mstroet/water_benchmarking/diagnostics
set -euo pipefail
gmx_mpi grompp -f rfinf.mdp -c start.gro -p spc.top -o rfinf.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm rfinf -ntomp 8 -nb gpu
printf "System\n" | gmx_mpi dipoles -f rfinf.xtc -s rfinf.tpr -o rfinf.M.xvg -eps rfinf.eps.xvg -epsilonRF 0 > rfinf.dip.log 2>&1
echo "DIAGNOSTICS3 COMPLETE"
