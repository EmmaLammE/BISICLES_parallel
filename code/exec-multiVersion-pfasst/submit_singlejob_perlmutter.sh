#!/bin/bash
#SBATCH --job-name=t64s64
#SBATCH --time=10:00:00
#SBATCH --qos=regular
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4096
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_mmp_t64s64.out

# unset unnecessary mem flag
echo "Before unset:"
env | grep -E 'SLURM_MEM' | sort
unset SLURM_MEM_PER_NODE
unset SLURM_MEM_PER_GPU
echo "After unset:"
env | grep -E 'SLURM_MEM' | sort

start_time=$(date +%s)
date +"Job started at: %Y-%m-%d %H:%M:%S"

# srun -n 1 ./driver2d.Linux.64.mpicxx.mpifort.DEBUG.MPI.ex ../data/inputs-noPython-pfasst.isomip24.melt4.l1l2.prlim.l0.Chombo.A2.2e-17.constfriction.sg4.a0.3-noPython-pfasst
srun -n 4096 ./driver2d.Linux.64.CC.mpifort.DEBUG.MPI.ex ../data/inputs-t64s64


# ---- resource report ----
sleep 60   # give Slurm a moment to record the finished srun step
echo "Nodes used: $SLURM_JOB_NUM_NODES -> $SLURM_JOB_NODELIST"
sacct -j $SLURM_JOB_ID --format=JobID%15,NNodes,NCPUS,NodeList%25,Elapsed,TotalCPU,MaxRSS,MaxRSSNode%12,State

# Estimated charge (pick the line for your cluster)
NCPU=$(sacct -n -X -j $SLURM_JOB_ID --format=NCPUS | tr -d ' ')
awk -v s=$SECONDS -v n=$SLURM_JOB_NUM_NODES 'BEGIN{printf "Perlmutter est: %.2f node-hours\n", s/3600*n}'
#awk -v s=$SECONDS -v c=$NCPU 'BEGIN{printf "Midway est: %.2f SUs\n", s/3600*c}'

end_time=$(date +%s)
date +"Job ended at: %Y-%m-%d %H:%M:%S"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

echo "Job took $hours hours, $minutes minutes, and $seconds seconds to complete."
