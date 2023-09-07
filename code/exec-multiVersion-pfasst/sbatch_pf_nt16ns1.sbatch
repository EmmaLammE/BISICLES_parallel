#!/bin/bash
#SBATCH --job-name=t16s1
#SBATCH --qos=regular
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=475G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_nt16ns1.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_nt16ns1.%j.err


srun ./driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
