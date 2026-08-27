#!/bin/bash
#SBATCH --account=m72
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --job-name=w_cpu
#SBATCH --output=%x-%j.out
module load gromacs/2024.3-mixed
export OMP_NUM_THREADS=4
cd /scratch/m72/mstroet/water_benchmarking/diagnostics
set -euo pipefail
# The same 2048 SPC box, reaction field at eps_rf = 61 and at infinity, on CPU only:
# does the finite-eps_rf dependence of the dipole fluctuation that GROMOS shows
# appear in GROMACS when the GPU kernel is out of the picture?
for tag in cpu61 cpuinf; do
  gmx_mpi grompp -f $tag.mdp -c start.gro -p spc.top -o $tag.tpr -maxwarn 3
  srun -N 1 -n 8 -c 4 gmx_mpi mdrun -deffnm $tag -ntomp 4 -nb cpu
  printf "System\n" | gmx_mpi dipoles -f $tag.xtc -s $tag.tpr -o $tag.M.xvg -eps $tag.eps.xvg -epsilonRF 61 > $tag.dip.log 2>&1 || true
done
echo "CPU DIAGNOSTICS COMPLETE"
