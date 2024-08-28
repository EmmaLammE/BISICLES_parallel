#!/bin/bash
#SBATCH --job-name=con3_n280
#SBATCH --time=1-00:00:00
#SBATCH --mem=40GB
#SBATCH -p serc
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH -o ./../../../sbatch_outputs/out_converge_nt280_lvl3.%j.out
#SBATCH -e ./../../../sbatch_outputs/err_converge_nt280_lvl3.%j.err

# below you run/call your code, load modules, python, Matlab, R, etc.
# and do any other scripting you want
# lines that begin with #SBATCH are directives (requests) to the scheduler-SLURM module load python/3.6.1
export LD_LIBRARY_PATH=/home/groups/jsuckale/liuwj/Toolkits/hdf5-1.12.2/lib:$LD_LIBRARY_PATH
srun -n 1 ./../../../driver2d.Linux.64.mpicxx.gfortran.DEBUG.MPI.ex inputs.test
cp ref_amr1_regrid10_ns1_nt280_tfin50_computeVel.000000.2d.hdf5 ref_start.hdf5
cp ref_amr1_regrid10_ns1_nt280_tfin50_computeVel.000280.2d.hdf5 ref_end.hdf5
rm ref_amr1_regrid10_ns1_nt280_tfin50_computeVel.*