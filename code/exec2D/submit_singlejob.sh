#!/bin/bash
#SBATCH --job-name=mmp
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=106GB
#SBATCH -p serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -o ./sbatch_outputs/out_pf_mmp_t1s1.out
#SBATCH -e ./sbatch_outputs/err_pf_mmp_t1s1.err


start_time=$(date +%s)
date +"Job started at: %Y-%m-%d %H:%M:%S"

srun -n 1 ./driver2d.Linux.64.mpicxx.mpifort.DEBUG.OPTHIGH.MPI.ex inputs.mismip


end_time=$(date +%s)
date +"Job ended at: %Y-%m-%d %H:%M:%S"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "Job took $hours hours, $minutes minutes, and $seconds seconds to complete."

