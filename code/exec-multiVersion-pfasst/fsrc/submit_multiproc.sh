#!/bin/bash

# Define arrays with number of processors for time and space dimensions
processors_time=(1 2)
processors_space=(1 2)

for nt in "${processors_time[@]}"; do
    for ns in "${processors_space[@]}"; do
        # total number of processors 
        np=$((nt * ns))

        # input file based on the number of processors in time and space
        input_file="./inputs_multiruns/inputs_t${nt}s${ns}.test"

        # Adjust job name and file outputs to include the number of processors
        job_name="t${nt}s${ns}"
        output_file="./sbatch_outputs/out_pf_amr3_t${nt}s${ns}.%j.out"
        error_file="./sbatch_outputs/err_pf_amr3_t${nt}s${ns}.%j.err"

        # temporary batch script for each job
        batch_script="temp_job_${nt}_${ns}.sh"
        cat > $batch_script << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --qos=regular
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$np
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --constraint=cpu
#SBATCH -o $output_file
#SBATCH -e $error_file

srun -n $np ./driver2d.Linux.64.CC.ftn.DEBUG.MPI.ex $input_file
EOF

        # Submit
        sbatch $batch_script

        # remove the temporary script
        rm $batch_script

    done
done