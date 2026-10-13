#!/bin/bash
#SBATCH --job-name=mmp
#SBATCH --time=3-00:00:00
#SBATCH -p serc
#SBATCH --mem-per-cpu=16GB
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH -o ./sbatch_outputs/out_pf_mmp_t1s1.out
#SBATCH -e ./sbatch_outputs/err_pf_mmp_t1s1.err

# unset unnecessary mem flag
echo "Before unset:"
env | grep -E 'SLURM_MEM' | sort
unset SLURM_MEM_PER_NODE
unset SLURM_MEM_PER_GPU
echo "After unset:"
env | grep -E 'SLURM_MEM' | sort

start_time=$(date +%s)
date +"Job started at: %Y-%m-%d %H:%M:%S"

srun -n 1 ./driver2d.Linux.64.mpicxx.mpifort.DEBUG.MPI.ex ../data/inputs-noPython-pfasst.isomip24.melt4.l1l2.prlim.l0.Chombo.A2.2e-17.constfriction.sg4.a0.3-noPython-pfasst

end_time=$(date +%s)
date +"Job ended at: %Y-%m-%d %H:%M:%S"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "Job took $hours hours, $minutes minutes, and $seconds seconds to complete."
