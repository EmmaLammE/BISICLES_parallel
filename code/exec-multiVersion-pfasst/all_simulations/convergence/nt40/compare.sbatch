#!/bin/bash
#SBATCH --job-name=com_nt40
#SBATCH --time=00:05:00
#SBATCH --mem=40GB
#SBATCH -p serc
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH -o ./../../../sbatch_outputs/out_convergeCompare_nt40.%j.out
#SBATCH -e ./../../../sbatch_outputs/err_convergeCompare_nt40.%j.err

# below you run/call your code, load modules, python, Matlab, R, etc.
# and do any other scripting you want
# lines that begin with #SBATCH are directives (requests) to the scheduler-SLURM module load python/3.6.1
export LD_LIBRARY_PATH+=/home/groups/jsuckale/liuwj/Toolkits/hdf5-1.12.2/lib:$LD_LIBRARY_PATH
srun -n 1 ./../../../ChomboCompare/compare2d.Linux.64.mpicxx.gfortran.DEBUG.MPI.ex inputs.compare.pf