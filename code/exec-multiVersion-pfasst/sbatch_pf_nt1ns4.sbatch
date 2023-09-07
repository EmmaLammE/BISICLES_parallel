#!/bin/bash
#SBATCH --job-name=t1s4
#SBATCH --qos=regular
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=475G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_nt1ns4.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_nt1ns4.%j.err


srun ./driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
