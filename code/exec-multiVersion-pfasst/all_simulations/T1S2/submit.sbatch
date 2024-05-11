#!/bin/bash
#SBATCH --job-name=t1s2_32
#SBATCH --qos=overrun
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --constraint=cpu
#SBATCH -o out_pf_nt1ns2_grid32_nt128.%j.out
#SBATCH -e err_pf_nt1ns2_grid32_nt128.%j.err


srun -n 2 ../.././driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
