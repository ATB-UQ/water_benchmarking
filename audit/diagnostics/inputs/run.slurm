#!/bin/bash
#SBATCH --account=m72-gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=02:00:00
#SBATCH --job-name=w_diag
#SBATCH --output=%x-%j.out
module load gromacs-amd-gfx90a/2026.1-mixed
export OMP_NUM_THREADS=8
cd /scratch/m72/mstroet/water_benchmarking/diagnostics
set -euo pipefail
for tag in vrescale pme; do
  gmx_mpi grompp -f $tag.mdp -c start.gro -p spc.top -o $tag.tpr -maxwarn 3
  srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest \
    gmx_mpi mdrun -deffnm $tag -ntomp 8 -nb gpu
done
echo "DIAGNOSTICS COMPLETE"
