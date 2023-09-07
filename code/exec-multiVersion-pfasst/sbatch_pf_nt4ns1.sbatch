#!/bin/bash
#SBATCH --job-name=t4s1
#SBATCH --qos=preempt
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=475G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_nt4ns1.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_nt4ns1.%j.err


srun ./driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex inputs.test
