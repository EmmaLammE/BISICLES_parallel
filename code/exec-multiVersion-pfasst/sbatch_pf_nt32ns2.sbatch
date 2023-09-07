#!/bin/bash
#SBATCH --job-name=t32s2
#SBATCH --qos=overrun
#SBATCH --time=00:10:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=475G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_nt32ns2.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_nt32ns2.%j.err


srun ./driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
