#!/bin/bash
#SBATCH --job-name=t2s16
#SBATCH --qos=regular
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=450G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_nt2ns16.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_nt2ns16.%j.err


srun ./driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
