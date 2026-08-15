#!/bin/bash

# Count the number of running jobs with the job name "minimize"
running_jobs=$(squeue -u $USER --name=minimize --state=R | wc -l)

# Count the number of pending jobs with the job name "minimize"
pending_jobs=$(squeue -u $USER --name=minimize --state=PD | wc -l)

# Adjust the counts since squeue includes a header line
running_jobs=$((running_jobs - 1))
pending_jobs=$((pending_jobs - 1))

# Output the results
echo "Number of running jobs with name 'minimize': $running_jobs"
echo "Number of pending jobs with name 'minimize': $pending_jobs"

