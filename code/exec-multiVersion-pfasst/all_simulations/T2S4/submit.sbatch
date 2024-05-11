#!/bin/bash
#SBATCH --job-name=t2s4
#SBATCH --qos=overrun
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --constraint=cpu
#SBATCH -o out_pf_nt2ns4.%j.out
#SBATCH -e err_pf_nt2ns4.%j.err


srun -n 8 ../.././driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
