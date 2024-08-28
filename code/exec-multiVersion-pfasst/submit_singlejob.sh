#!/bin/bash
#SBATCH --job-name=t4s4
#SBATCH --qos=overrun
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_amr2_nt2000_tfin20_t4s4.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_amr2_nt2000_tfin20_t4s4.%j.err

start_time=$(date +%s)
date +"Job started at: %Y-%m-%d %H:%M:%S"

# Run the main application
srun -n 16 ./driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex ./inputs_multiruns/inputs_nt2000_tfin20_t4s4.test

end_time=$(date +%s)
date +"Job ended at: %Y-%m-%d %H:%M:%S"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "Job took $hours hours, $minutes minutes, and $seconds seconds to complete."

