#!/bin/bash
#SBATCH --job-name=a2nt4000
#SBATCH --partition=serc
#SBATCH --time=4-15:00:00
#SBATCH --ntasks-per-node=1       # Example configuration
#SBATCH --cpus-per-task=10          # This configuration uses exactly 128 CPUs per node
#SBATCH --mem-per-cpu=10G                   # Assuming this memory is sufficient per node
#SBATCH -o ./sbatch_outputs/out_pf_amr2_nt4000_tfin20.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_amr2_nt4000_tfin20.%j.err
# export SLURM_CPUS_PER_TASK=10

export SLURM_CPU_BIND=none

start_time=$(date +%s)
date +"Job started at: %Y-%m-%d %H:%M:%S"

srun ./driver2d.Linux.64.mpicxx.mpifort.DEBUG.OPTHIGH.MPI.ex ./inputs_multiruns/inputs_nt4000_tfin20_t1s1.test

end_time=$(date +%s)
date +"Job ended at: %Y-%m-%d %H:%M:%S"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "Job took $hours hours, $minutes minutes, and $seconds seconds to complete."

