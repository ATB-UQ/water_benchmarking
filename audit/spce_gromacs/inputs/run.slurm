#!/bin/bash
#SBATCH --account=m72-gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=06:00:00
#SBATCH --job-name=w_spce
#SBATCH --output=%x-%j.out

module load gromacs-amd-gfx90a/2026.1-mixed

# mdrun's -ntomp must agree with what SLURM allocated, or it oversubscribes.
export OMP_NUM_THREADS=8
cd /scratch/m72/mstroet/water_benchmarking/spce
set -euo pipefail

gmx_mpi grompp -f emin.mdp -c water_2048.gro -p spce.top -o emin.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm emin -ntomp 8

gmx_mpi grompp -f eq1.mdp -c emin.gro -p spce.top -o eq1.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm eq1 -ntomp 8 -nb gpu

gmx_mpi grompp -f eq2.mdp -c eq1.gro -p spce.top -o eq2.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm eq2 -ntomp 8 -nb gpu

gmx_mpi grompp -f eq3.mdp -c eq2.gro -p spce.top -o eq3.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm eq3 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_01.mdp -c eq3.gro -p spce.top -o md_01.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_01 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_02.mdp -c md_01.gro -p spce.top -o md_02.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_02 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_03.mdp -c md_02.gro -p spce.top -o md_03.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_03 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_04.mdp -c md_03.gro -p spce.top -o md_04.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_04 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_05.mdp -c md_04.gro -p spce.top -o md_05.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_05 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_06.mdp -c md_05.gro -p spce.top -o md_06.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_06 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_07.mdp -c md_06.gro -p spce.top -o md_07.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_07 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_08.mdp -c md_07.gro -p spce.top -o md_08.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_08 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_09.mdp -c md_08.gro -p spce.top -o md_09.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_09 -ntomp 8 -nb gpu

gmx_mpi grompp -f md_10.mdp -c md_09.gro -p spce.top -o md_10.tpr -maxwarn 3
srun -l -u -N 1 -n 1 -c 8 --gpus-per-node=1 --gpus-per-task=1 --gpu-bind=closest gmx_mpi mdrun -deffnm md_10 -ntomp 8 -nb gpu

echo "LADDER COMPLETE"