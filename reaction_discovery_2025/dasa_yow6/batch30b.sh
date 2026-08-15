#!/usr/bin/bash

for i in {1..30}; do
    dir="secondary$i"
    cd "$dir"
    
    # Submit the job and capture the job ID
    job_output=$(sbatch submit.sh)
    job_id=$(echo "$job_output" | awk '{print $4}')
    echo "Submitted job $job_id in directory $dir"

    # Check the job status until it's running ("R")
    while true; do
        job_status=$(squeue -j "$job_id" -o "%T" -h)

        if [[ "$job_status" == "RUNNING" ]]; then
            echo "Job $job_id is now running."
            break
        else
            echo "Job $job_id is currently $job_status. Waiting for it to start running..."
            sleep 5  # Check every 5 seconds
        fi
    done
    
    cd ..
done

