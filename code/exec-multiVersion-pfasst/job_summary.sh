#!/bin/bash
job_id=$1

# Use sacct to get detailed job information
sacct -j $job_id --format=JobID,JobName,State,Elapsed,TotalCPU,NCPUS,MaxRSS,AveRSS,ReqMem,CPUTime,SystemCPU > job_summary.txt

# Extract details for the main job step (e.g., job step .0)
line=$(grep "\.0 " job_summary.txt)
echo $line
elapsed=$(echo $line | awk '{print $4}')
total_cpu=$(echo $line | awk '{print $5}')
ncpus=$(echo $line | awk '{print $6}')
averss_kb=$(echo $line | awk '{print $8}' | sed 's/K//')

# Convert elapsed time to seconds
IFS=: read h m s <<< "$elapsed"
elapsed_sec=$((10#$h*3600 + 10#$m*60 + 10#$s))

# Convert total CPU time to seconds
IFS=: read h m s <<< "$total_cpu"
total_cpu_sec=$((10#$h*3600 + 10#$m*60 + 10#$s))

# Convert Avg RSS to MB and GB
# echo $averss_kb
# echo $($averss_kb/1024)
if [[ "$averss_kb" =~ ^[0-9]+$ ]]; then
    # Convert Avg RSS to MB and GB
    averss_mb=$(awk "BEGIN {printf \"%.2f\", $averss_kb / 1024}")
    averss_gb=$(awk "BEGIN {printf \"%.2f\", $averss_kb / 1048576}")
else
    averss_mb="N/A"
    averss_gb="N/A"
fi

# Calculate CPU efficiency, ensuring no division by zero
if [ "$elapsed_sec" -gt 0 ] && [ "$ncpus" -gt 0 ]; then
    cpu_efficiency=$(awk "BEGIN {printf \"%.2f\", ($total_cpu_sec / ($elapsed_sec * $ncpus)) * 100}")
else
    cpu_efficiency="N/A"
fi

echo "Job Summary for $job_id"
echo "Elapsed Time: $elapsed"
echo "Total CPU Time: $total_cpu"
echo "Number of CPUs: $ncpus"
echo "CPU Efficiency: $cpu_efficiency%"
echo "Avg RSS: $averss_kb KB ($averss_mb MB, $averss_gb GB)"
