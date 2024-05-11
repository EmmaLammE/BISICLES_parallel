#!/bin/bash
#SBATCH --job-name=t8s1
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --constraint=cpu
#SBATCH -o out_pf_nt8ns1.%j.out
#SBATCH -e err_pf_nt8ns1.%j.err


srun -n 8 ../.././driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
