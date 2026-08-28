#!/bin/bash
#SBATCH --account=m72-gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=04:00:00
#SBATCH --job-name=w_o3set
#SBATCH --output=%x-%j.out
module load gromacs-amd-gfx90a/2026.1-mixed
export OMP_NUM_THREADS=8
cd /scratch/m72/mstroet/water_benchmarking/opc3_settings
set -euo pipefail

# One 100 ps NPT re-equilibration at the new cutoff and timestep, from the box the
# 1.8 nm / 1 fs production run finished on.  The full ladder is not repeated -- the
# box is already liquid water at 298 K -- but the equilibrium density does depend on
# the cutoff, so going straight to production would put that relaxation inside the
# averages.
gmx_mpi grompp -f o3eq.mdp -c start.gro -p opc3.top -o o3eq.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest \
  gmx_mpi mdrun -deffnm o3eq -ntomp 8 -nb gpu

# The three conditions all branch from that same re-equilibrated box, so the only
# difference between them is the electrostatics.
for tag in o3rf78 o3rf61 o3pme; do
  gmx_mpi grompp -f $tag.mdp -c o3eq.gro -p opc3.top -o $tag.tpr -maxwarn 3
  # -pme gpu is fatal under a reaction field, so only the PME run offloads it.
  extra=""
  if [ "$tag" = "o3pme" ]; then extra="-pme gpu"; fi
  srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest \
    gmx_mpi mdrun -deffnm $tag -ntomp 8 -nb gpu $extra
  grep -A2 "Performance:" $tag.log | tail -2
done
echo "OPC3 SETTINGS SWEEP COMPLETE"
