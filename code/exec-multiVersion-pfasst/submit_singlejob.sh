#!/bin/bash
#SBATCH --job-name=t1s1
#SBATCH --qos=regular
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_amr2_nt1200_tfin20_t1s1.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_amr2_nt1200_tfin20_t1s1.%j.err


srun ./driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex inputs.test
