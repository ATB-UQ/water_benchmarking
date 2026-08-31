#!/bin/bash
#SBATCH --account=m72-gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=05:00:00
#SBATCH --job-name=w_diag2
#SBATCH --output=%x-%j.out
module load gromacs-amd-gfx90a/2026.1-mixed
export OMP_NUM_THREADS=8
cd /scratch/m72/mstroet/water_benchmarking/diagnostics
set -euo pipefail
RUN="srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -ntomp 8 -nb gpu"
# Same 2048 box, shorter cutoffs: Rc/L = 0.35 and 0.23
for tag in rc14 rc09; do
  gmx_mpi grompp -f $tag.mdp -c start.gro -p spc.top -o $tag.tpr -maxwarn 3
  $RUN -deffnm $tag
done
# 8x box (2x2x2 replication of the equilibrated 2048 box), 1.8 nm cutoff: Rc/L = 0.23
gmx_mpi genconf -f start.gro -nbox 2 2 2 -o big.gro
gmx_mpi grompp -f big18.mdp -c big.gro -p spc_big.top -o big18.tpr -maxwarn 3
$RUN -deffnm big18
echo "DIAGNOSTICS2 COMPLETE"
