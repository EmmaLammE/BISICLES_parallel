#!/bin/bash
#SBATCH --job-name=sc_t2s1
#SBATCH --time=02:00:00
#SBATCH --mem=40GB
#SBATCH -p serc
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=2
#SBATCH -o out_scale_t2s1.%j.out
#SBATCH -e err_scale_t2s1.%j.err

# below you run/call your code, load modules, python, Matlab, R, etc.
# and do any other scripting you want
# lines that begin with #SBATCH are directives (requests) to the scheduler-SLURM module load python/3.6.1
export LD_LIBRARY_PATH=/home/groups/jsuckale/liuwj/Toolkits/hdf5-1.12.2/lib:$LD_LIBRARY_PATH
srun --cpus-per-task=1 -n 2 ./../../../driver2d.Linux.64.mpicxx.gfortran.DEBUG.MPI.ex inputs.test
cp ref_amr1_regrid10_ns2_nt128_tfin50_computeVel.000000.2d.hdf5 ref_start.hdf5
cp ref_amr1_regrid10_ns2_nt128_tfin50_computeVel.000128.2d.hdf5 ref_end.hdf5
rm ref_amr1_regrid10_ns2_nt128_tfin50_computeVel.*