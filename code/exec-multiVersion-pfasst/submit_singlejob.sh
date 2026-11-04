#!/bin/bash
#SBATCH --job-name=t2s2
#SBATCH --qos=regular
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --constraint=cpu
#SBATCH -o ./sbatch_outputs/out_pf_twisty_test.%j.out
#SBATCH -e ./sbatch_outputs/err_pf_twisty_test.%j.err


#srun ./driver2d.Linux.64.g++.gfortran.DEBUG.MPI.ex ../data/inputs-noPython-pfasst.isomip24.melt4.l1l2.prlim.l0.Chombo.A2.2e-17.constfriction.sg4.a0-1.3-noPython-pfasst

srun ./driver2d.Linux.64.g++.gfortran.DEBUG.MPI.ex inputs_multirun/inputs_t2s2.test


