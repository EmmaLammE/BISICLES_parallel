#!/bin/bash
#SBATCH --job-name=t2s2
#SBATCH --qos=regular
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_amr2_nt4_tfin1_t2s2.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_amr2_nt4_tfin1_t2s2.%j.err


srun -n 4 ./driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex inputs.test
