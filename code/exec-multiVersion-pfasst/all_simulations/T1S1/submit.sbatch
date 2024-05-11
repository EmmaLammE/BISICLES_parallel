#!/bin/bash
#SBATCH --job-name=t1s1
#SBATCH --qos=overrun
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --constraint=cpu
#SBATCH -o out_pf_nt1ns1.%j.out
#SBATCH -e err_pf_nt1ns1.%j.err


srun -n 1 ../.././driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
