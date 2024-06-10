#!/bin/bash
#SBATCH --job-name=t8s32
#SBATCH --qos=regular
#SBATCH --time=03:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_amr2_nt256_tfin20_t8s32.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_amr2_nt256_tfin20_t8s32.%j.err


srun -n 256 ./driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex ./inputs_multiruns/inputs_nt256_tfin20_t8s32.test