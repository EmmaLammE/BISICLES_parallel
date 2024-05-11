#!/bin/bash
#SBATCH --job-name=t2s1
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --constraint=cpu
#SBATCH -o out_pf_nt1ns4.%j.out
#SBATCH -e err_pf_nt1ns4.%j.err


srun -n 4 ../.././driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
