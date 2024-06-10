#!/bin/bash

# Define arrays with number of time steps
Nt=(140 150 160 170 180 190 200)

for nt in "${Nt[@]}"; do
    # input file based on the number of time steps
    input_file="./inputs_multiruns/inputs_nt${nt}_tfin20_t1s1.test"

    # Adjust job name and file outputs to include the number of processors
    job_name="nt${nt}"
    output_file="./sbatch_outputs/out_pf_amr2_nt${nt}_tfin20_t1s1.%j.out"
    error_file="./sbatch_outputs/err_pf_amr2_nt${nt}_tfin20_t1s1.%j.err"

    # temporary batch script for each job
    batch_script="temp_job_nt${nt}.sh"
    cat > $batch_script << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --qos=regular
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --constraint=cpu
#SBATCH -o $output_file
#SBATCH -e $error_file
srun ./driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.ex $input_file
EOF

    # Submit job and capture the job ID
    job_id=$(sbatch $batch_script | cut -d " " -f 4)
    echo "Submitted job $job_id"

    # Create and submit a dependent cleanup job
    cleanup_script="cleanup_nt${nt}.sh"
    cat > $cleanup_script << EOF
#!/bin/bash
#SBATCH --job-name=clean_nt${nt}
#SBATCH --dependency=afterok:$job_id
#SBATCH --qos=regular
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --constraint=cpu
# keep only the first and last data files

prefix="ref_amr2_nt\${nt}_tfin20"
declare -a indices
for file in \${prefix}.*.hdf5; do
    index=\$(echo "\$file" | sed "s/\${prefix}\.\(.*\)\.hdf5/\1/")
    indices+=("\$index")
done
max_index=\$(printf "%s\n" "\${indices[@]}" | sort -n | tail -1)
for file in \${prefix}.*.hdf5; do
    index=\$(echo "\$file" | sed "s/\${prefix}\.\(.*\)\.hdf5/\1/")
    if [[ "\$index" != "000000" && "\$index" != "\$max_index" ]]; then
        rm "\$file"
    fi
done
EOF
    # Submit the cleanup job
    sbatch $cleanup_script

    # Remove the temporary scripts
    rm $batch_script $cleanup_script
done